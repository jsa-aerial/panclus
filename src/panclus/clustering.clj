(ns panclus.clustering
  [:require
   [clojure.data.csv :as csv]
   [clojure.string :as cljstr]

   [aerial.fs :as fs]

   [aerial.utils.string :as str]
   [aerial.utils.io :refer [letio] :as io]
   [aerial.utils.coll :refer [vfold] :as coll]
   [aerial.utils.math :as m]
   [aerial.utils.math.probs-stats :as p]
   [aerial.utils.math.infoth :as it]
   [aerial.utils.math.combinatorics :as cmb]

   [aerial.bio.utils.files :as bufiles]
   [aerial.bio.utils.aligners :as aln]
   [aerial.bio.utils.seqxform :as sqx]

   [panclus.params :as pams]
   ])




(def spsbase (pams/get-params :spsbase))

;;; Separate out the strains and aggregate the annotations per strain
;;; from original db querys

;;; (separate-strains "orig-db-genes-locus-tags-pg320.csv" "PG/PG320")
(defn separate-strains [orig-db-csv strain-dir]
  (letio [id-nm-map {"17" "locus_tag", "23" "protein_id",
                     "32" "gene", "47" "old_locus_tag"}
          tags ["gene" "old_locus_tag" "locus_tag" "protein_id"]
          cols (coll/concatv ["ACCEN" "Start" "End" "Strand"] tags)
          strainbase (fs/join spsbase strain-dir)
          dbcsv (fs/join strainbase orig-db-csv)
          lines (io/read-lines dbcsv)
          strains (group-by first (mapv #(-> % csv/read-csv first) lines))
          NCs (keys strains)]
    (doseq [NC NCs]
      (let [recs (strains NC)
            loci (group-by #(->> % rest (take 3) (cljstr/join "/")) recs)
            recs (map (fn[[_ v]]
                        (let [[nm s e st :as evec] (->> v first (take 4))
                              tag-map (into {} (map (fn[[_ _ _ _ nm id]]
                                                      [(id-nm-map id) nm]) v))
                              tags (mapv #(tag-map % "NA") tags)]
                          (coll/concatv evec tags)))
                      loci)
            ot (io/open-file (fs/join strainbase "CSV" (str NC ".csv")) :out)]
        (csv/write-csv ot (cons cols recs))
        (.close ot)))))

;;; Generate fastas from above generated CSVs
;;; (future  (gen-strain-fastas "PG/PG320"))
(defn gen-strain-fastas [strain-dir]
  (let [base spsbase
        strainbase (fs/join base strain-dir)
        tags ["gene" "old_locus_tag" "locus_tag" "protein_id"]
        csvs (fs/re-directory-files (fs/join strainbase "CSV") "*.csv")
        dbresolver (aerial.bio.utils.params/get-params [:genomes :tvoseq02])]
    (aerial.bio.utils.params/with-genome-db dbresolver
      (doseq [csv (sort csvs)]
        (letio [nc (->> csv fs/basename (str/split #"\.") first)
                fna (fs/join strainbase (str nc ".fna"))
                lines (io/read-lines csv)
                recs (rest (mapv #(-> % csv/read-csv first) lines))]
          (println fna)
          (io/with-out-writer fna
            (doseq [[nm s e st & tagvals] recs]
              (let [e (Integer. e)
                    cnt (inc (- (Integer. e) (Integer. s))) ;;_ (println cnt)
                    r (second (aerial.utils.math/div cnt 3))
                    e (cond (= r 0) e, (= r 1) (dec e), (= r 2) (inc e))
                    entry (aerial.bio.utils.files/make-entry nm s e st)
                    tagstr (cljstr/join "," (interleave tags tagvals))]
                (println (str ">" entry "," tagstr))
                (println (-> entry aerial.bio.utils.files/gen-name-seq second
                             aerial.bio.utils.seqxform/ntseq->aaseq))))))))))


(def pgfnas
  (let [nmfn #(-> % fs/basename (fs/replace-type ""))
        pg30fnas  (->> (fs/directory-files (fs/join spsbase "PG/PG30") "fna")
                       sort vec)
        pg320fnas (->> (fs/directory-files (fs/join spsbase "PG/PG320") "fna")
                      sort vec)]
    {:pg30 pg30fnas
     :pg320 pg320fnas}))

(def strains
  (let [nmfn #(-> % fs/basename (fs/replace-type ""))]
    {:pg30 (mapv nmfn (:pg30 pgfnas)) :pg320 (mapv nmfn (:pg320 pgfnas))}))



(defn ntseq2aaseq [s]
  (let [[_ r] (m/div (count s) 3)
        s (if (= 0 r) s (str/butlast r s))]
    (sqx/ntseq->aaseq s)))


(def ows {:nt 9 :aa 6})

(defn gen-sqs-dists
  "Generate the kmer distributions for each sq in 'sqs', where sq is
  the pair [entry ntsq] 'ows' is the optimal word size (from cre) and
  'aa' is boolean for whether should use ammino acid (convert first)"
  [sqs & {:keys [aa ows] :or {aa false ows 9}}]
  (vfold (fn[[n s]]
           (let [s (if aa (ntseq2aaseq s) s)])
           [n (count s) (p/probs ows s)])
         sqs))

(defn gen-strain-dists
  "Take a fasta file 'strain-fna', which contains the CDS/genes (as NT
  sqs) for one or more strains and produce a set of corresponding
  probability distributions for each CDS/gene using 'ows' kmer
  size. If 'aa' is true, first convert to ammino acid
  sequences. Returns a vector of triplets [ent len dist], one for each
  CDS/gene, where 'ent' is the entry defining the loci, 'len' is the
  length of the sq, and 'dist' is the probability distribution."
  [strain-fna & {:keys [aa ows] :or {aa false ows 9}}]
  (let [name-seq-pairs (bufiles/read-seqs strain-fna :info :both)]
    (gen-sqs-dists
     (sort-by second (fn[l r] (< (count l) (count r))) name-seq-pairs)
     :aa aa :ows ows)))

(defn gen-strain-group-dists
  "Take a set of fastas 'strain-fnas', and perform 'gen-strain-dists'
  on each of them. Useful for operating over a directory of such
  fastas."
  [strain-fnas & {:keys [aa ows] :or {aa false ows 9}}]
  (mapv #(gen-strain-dists % :ows ows :aa aa) strain-fnas))


(comment
  (def sqs-map (load-sqs))

  (def ref-dists-all
    (vfold (fn[[n s]] [n (count s) (p/probs (ows :nt) s)])
           (sqs-map :refsqs)))

  (def maxNCdists
    (filterv (fn[[n cnt P]]
               (= "NC_014498" (str/substring n 0 9)))
             ref-dists-all))
  (def ref-dists
    (filterv (fn[[n cnt P]]
               (not= "NC_014498" (str/substring n 0 9)))
             ref-dists-all))
  (def refdists-map
    (reduce (fn[M [n cnt P :as v]] (assoc M n v)) {} ref-dists))

  (def ref-len-dists-map
    (reduce (fn[S [n cnt P :as rec]]
              (assoc S cnt (conj (get S cnt []) rec)))
            {} ref-dists))
  (io/with-out-writer "/store/data/PanClus/Entropy/ref-len-dists-map.clj"
    (prn ref-len-dists-map))


  (let [base "/store/data/PanClus/Entropy"
        data-nms [[maxNCdists "maxNCdists.clj"]
                  [ref-dists "ref-dists.clj"]
                  [refdists-map "refdists-map.clj"]]]
    (doseq [[data nm] data-nms]
      (io/with-out-writer (fs/join base nm)
        (prn data))))

  (def fut
    (future
      (do (def pg30dists
            (vfold (fn[[n s]] [n (count s) (p/probs (ows :nt) s)])
                   (sqs-map :pg30sqs)))
          (def pg320dists
            (vfold (fn[[n s]] [n (count s) (p/probs (ows :nt) s)])
                   (sqs-map :pg320sqs))))))

  (def pg30dist-map
    (reduce (fn[M [n cnt P :as v]] (assoc M n v)) {} pg30dists))
  (def pg320dist-map
    (reduce (fn[M [n cnt P :as v]] (assoc M n v)) {} pg320dists))

  (let [base "/store/data/PanClus/Entropy"
        data-nms [[pg30dists "pg30dists.clj"]
                  [pg320dists "pg320dists.clj"]
                  [pg30dist-map "pg30dists-map.clj"]
                  [pg320dist-map "pg320dists-map.clj"]]]
    (doseq [[data nm] data-nms]
      (io/with-out-writer (fs/join base nm)
        (prn data))))

  )



(def ref-len-dists-map
  (->>  "/store/data/PanClus/Entropy/ref-len-dists-map.clj"
        slurp read-string))
(def maxNCdists
  (->> "/store/data/PanClus/Entropy/maxNCdists.clj"
       slurp read-string))






(defn opt-wz [name-seq-pairs &
              {:keys [limit cnt crecut alpha]
               :or {limit 14 cnt 5 crecut 0.10}}]
  (let [cres (for [[nm sq] name-seq-pairs]
               (vfold (fn[l] [l (it/CREl l sq :alpha alpha)
                             (count sq) nm])
                      (range 3 (inc limit))))
        wz (reduce (fn [res v]
                     (/ (+ res
                           (some #(when (< (second %) crecut) (first %)) v))
                        2.0))
                   0.0 cres)]
    [(Math/round wz) cres]))

#_(filterv (fn[[n cnt Q]]
             (if (< mean-std pcnt mean+std)
               (< (Math/abs (- pcnt cnt)) len-max)
               (>= (-> (- 1 (/ (Math/abs (- pcnt cnt))
                               (double (max pcnt cnt))))
                       (* 10.0) Math/round (/ 10.0))
                   %-max)))
           (vals distmap))
(defn ld-map-mean [ld-map]
  (->> ld-map
     (reduce (fn[[n d] [k v]]
               (let [cnt (count v)] [(+ n (* cnt k)) (+ d cnt)]))
             [0.0 0.0])
     (apply /)))

(defn cluster-strains
  [dists distmap maxNCdists ref-len-dists-map
   & {:keys [len-max %-max jsdctpt] :or {len-max 5 %-max 0.99 jsdctpt 0.40}}]

  (let [mean (p/mean (map second dists))
        std  (p/std-deviation (map second dists))
        mean-std (long (- mean std))
        mean+std (long (+ mean std 7))]
    (loop [;;distmap distmap
           maxNCdists maxNCdists
           clus []]
      (if (empty? maxNCdists)
        clus
        (let [[n1 pcnt P] (first maxNCdists)
              ds (if (< mean-std pcnt mean+std)
                   (->> (range (- pcnt len-max -1)
                               (+ pcnt len-max))
                        (mapcat ref-len-dists-map))
                   (->> (range (-> pcnt (* 0.99) Math/round)
                               (-> pcnt (* 1.01) Math/round))
                        (mapcat ref-len-dists-map)))
              ;;_ (println :ds-cnt (count ds))
              clu (coll/takev-while
                   #(-> % second (< jsdctpt))
                   (sort-by
                    second (vfold (fn[[n cnt Q]]
                                    [n (it/jensen-shannon P Q)])
                                  ds)))
              ;;_ (println :clu-cnt (count clu))
              ;;newdists (apply dissoc distmap (mapv first clu))
              ;;_ (println :newdist-cnt (count newdists))
              clu (conj clu [n1 pcnt 0.0])]
          (recur #_newdists
                 (rest maxNCdists)
                 (conj clus clu)))))))

(defn make-hybrids
  "Transform the clusters in 'clusters' obtained from distribution map
  'distmap' into dist recs of form: [\"hybrid-dict-x\", cnt,
  hybrid-dict] where 'cnt' is the mean seq lens in cluster x and
  'hybrid-dict' is the minimum entropy derived hybrid dictionary /
  distribution of all distributions in cluster x"
  [clusters distmap & {:keys [wz] :or {wz 6}}]
  (map
   (fn[clu i]
     (let [recs (->> clu (map first) (map distmap))
           dists (mapv last recs)
           cnt (long (p/mean (mapv second recs)))]
       [(format "hbrid-dict-%s" (inc i)) cnt (it/hybrid-dictionary 6 dists)]))
   clusters (range)))

(defn cluster-leftovers
  [distmap dists all-matched ref-len-dists-map
   & {:keys [len-max %-max jsdctpt] :or {len-max 5 %-max 0.99 jsdctpt 0.40}}]
  (let [distmap (reduce (fn[M k] (dissoc M k)) distmap all-matched)
        ref-len-dists-map (reduce (fn[S [n cnt P :as rec]]
                               (assoc S cnt (conj (get S cnt []) rec)))
                             {} (vals distmap))
        mean (p/mean (map second dists))
        std  (p/std-deviation (map second dists))
        mean-std (long (- mean std))
        mean+std (long (+ mean std 7))]
    (loop [distmap distmap
           clus []]
      (if (empty? distmap)
        clus
        (let [[n1 pcnt P] (->> distmap first second)
              distmap (dissoc distmap n1)
              ds (if (< mean-std pcnt mean+std)
                   (->> (range (- pcnt len-max -1)
                               (+ pcnt len-max))
                        (mapcat ref-len-dists-map))
                   (->> (range (-> pcnt (* 0.99) Math/round)
                               (-> pcnt (* 1.01) Math/round))
                        (mapcat ref-len-dists-map)))
              ;;_ (println :ds-cnt (count ds))
              clu (coll/takev-while
                   #(-> % second (< jsdctpt))
                   (sort-by
                    second (vfold (fn[[n cnt Q]]
                                    [n (it/jensen-shannon P Q)])
                                  ds)))
              newdists (apply dissoc distmap (mapv first clu))
              ;;_ (println (count newdists))
              clu (conj clu [n1 pcnt 0.0])]
          (recur newdists
                 (conj clus clu)))))))




(def ref-clusters-040
  (time (cluster-strains ref-dists refdists-map
                         maxNCdists ref-len-dists-map :jsdctpt 0.40)))
(def all-matched-040
  (->> ref-clusters-040 (mapcat (fn[clu](->> clu butlast (map first))))))
(io/with-out-writer
  "/store/data/PanClus/Entropy/ref-clusters-040.clj"
  (prn ref-clusters-040))

(def ref-clu-map040
  (reduce (fn[M clu]
            (assoc M (->> clu last first (str/split #",") first) clu))
          {} ref-clusters-040))
(io/with-out-writer
  "/store/data/PanClus/Entropy/ref-clu-map-040"
  (prn ref-clu-map040))


(def ref-clus-leftovers-040
  (time (cluster-leftovers
         refdists-map ref-dists all-matched-040 ref-len-dists-map)))

(def ref-clus-leftover-map040
  (reduce (fn[M clu]
            (assoc M (->> clu last first (str/split #",") first) clu))
          {} ref-clus-leftovers-040))
(io/with-out-writer
  "/store/data/PanClus/Entropy/ref-clu-leftover-map-040"
  (prn ref-clus-leftover-map040))

(def total-ref-clus-map040
  (merge ref-clu-map040 ref-clus-leftover-map040))
(io/with-out-writer
  "/store/data/PanClus/Entropy/total-ref-clus-map-040"
  (prn total-ref-clus-map040))

(map (fn[cnt]
       [cnt (->> ref-clus-leftovers-040 (filter #(= cnt (count %))) count)])
     (range 1 34))


(->> ref-clusters-040 (map count) (filter #(<= 19 % 42)) count)
(->> ref-clusters-040 (map count) (filter #(< % 19)) count)
(->> ref-clusters-040 (map count) (filter #(< % 19)) (p/mean))
(->> ref-clusters-040 (map count) (filter #(< % 19)))

(def ref-clu-setmap040
  (->> ref-clusters-040
       (map (fn[clu]
              (let [cluset (->> clu (map first) (map #(str/split #"," %))
                                (map first) (map bufiles/entry-parts)
                                (map first) set)
                    seed (->> clu last first (str/split #",") first)]
                [seed cluset])))
       (into {})))
(io/with-out-writer
  "/store/data/PanClus/Entropy/ref-clu-setmap-040"
  (prn ref-clu-setmap040))
(def foo (->> ref-clu-setmap040 (filter #(<= 19 (->> % second count) 19))))



(->> ref-clusters-040 (map count) (filter #(> % 21)))
(def big1
  (->> ref-clusters-040 (filter #(> (count %) 21))
       (drop-while #(< (count %) 100)) first))
(->> big1 first first (str/split #",") first bufiles/entry-parts second)

(->> big1 (drop 50) first first (str/split #",")
     first bufiles/gen-name-seq second
     aerial.bio.utils.seqxform/ntseq->aaseq)
(->> big1 first first (str/split #",")
     first bufiles/gen-name-seq second
     aerial.bio.utils.seqxform/ntseq->aaseq)
(->> big1 (drop 25) first first (str/split #",")
     first bufiles/gen-name-seq second
     aerial.bio.utils.seqxform/ntseq->aaseq)

(->> big1 (drop-while #(< (last %) 0.2)) first first
     (str/split #",") first bufiles/gen-name-seq second
     aerial.bio.utils.seqxform/ntseq->aaseq)

(reduce (fn[M [id _ _]]
          (let [id (->> id (str/split #",") first bufiles/entry-parts first)]
            (assoc M id (inc (get M id 0)))))
        {} big1)

(->> big1
     (keep (fn[[id _ _]]
             (let [[id se st] (->> id (str/split #",") first
                                   bufiles/entry-parts)]
               (when (= id "NC_017591") [id se st]))))
     (sort-by #(-> % second first)))


(defn aln-score-clu [clu & {:keys [full?]}]
  (let [nt2aa aerial.bio.utils.seqxform/ntseq->aaseq
        seed-ent (->> clu last first (str/split #",") first)
        seed-aasq (->> seed-ent bufiles/gen-name-seq
                       second aerial.bio.utils.seqxform/ntseq->aaseq)
        best (count seed-aasq)
        clucnt (count clu)
        score-vecs
        (->> clu butlast
             (map #(let [ent (->> % first (str/split #",") first)
                         jsd (last %)
                         aasq (->> ent bufiles/gen-name-seq
                                   second nt2aa)]
                     [ent jsd aasq]))
             (map (fn[[ent jsd aasq]]
                    [ent jsd
                     (aln/hirsch-align
                      seed-aasq aasq
                      :match 1 :mmatch -1 :gap -1 :score-only? true)])))
        scores (mapv last score-vecs)
        [maxaln minaln meanaln] (if (seq scores)
                                  (map #(apply % scores)
                                       [max min p/mean])
                                  [best best best])
        max% (double (/ maxaln best))
        min% (double (/ minaln best))
        mean% (double (/ meanaln best))]
    [seed-ent clucnt best maxaln minaln meanaln max% min% mean%
     (if full? score-vecs)]))

(aln-score-clu big1)
(def cluscores
  (future (->> ref-clusters-040 (coll/drop-until #(> (count %) 1))
               (mapv aln-score-clu))))
(io/with-out-writer
  "/store/data/PanClus/Entropy/aln-scores.clj"
  (prn @cluscores))

(def clualn<080
  (->> @cluscores
       (sort-by #(->> % butlast last))
       (take (->> @cluscores (map #(->> % butlast last))
                  sort (take-while #(< % 0.8)) count))))
(->> clualn<080 (map #(ref-clu-map040 (first %))) first)
(->> clualn<080 (map #(ref-clu-map040 (first %))) (drop 3) first)


(->> ref-clusters (filter #(> (count %) 42)) first)
(->> ref-clusters (filter #(> (count %) 42)) first butlast last)
(->> ref-clusters
     (map (fn[clu] (if (= 1 (count clu)) 0.0 (->> clu butlast last last))))
     (sort >) (take-while #(> % 0.45)) count)


(cluster-strains ref-dists refdists-map (take 5 (coll/dropv 1000 maxNCdists)))
(cluster-strains ref-dists refdists-map (take 5 (coll/dropv 2000 maxNCdists)))


(def size-dist
  (sort-by
   first
   (reduce (fn[S [n cnt P]] (assoc S cnt (inc (get S cnt 0)))) {} ref-dists)))

(p/mean (map second ref-dists))
281.6591802747837
(p/median (map second ref-dists))
238.0
(p/std-deviation (map second ref-dists))
221.6938718239762
(- 282 222)
60
(+ 282 222)
504
(->> maxNCdists (take-while #(< (second %) 60)) count)
120
(->> maxNCdists (take-while #(< (second %) 504)) count)

(def jsadists
  (conj (map refdists-map (butlast (map first (first *1))))
        (first (coll/dropv 1000 maxNCdists))))

(p/mean (map second (map refdists-map (map first (rest jsadists)))))

(def hybdist
  ["hybrid-dist-test" 200 (it/hybrid-dictionary 0 (map last jsadists))])

(count (last hybdist))
236
(map count (map last (map refdists-map (map first (rest jsadists)))))
(195 195 195 195 195 195 195 195 195 195 195 195 194)

(cluster-strains ref-dists refdists-map [hybdist])
