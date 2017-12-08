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




(def panclus-base (pams/get-params :panclus-base))

;;; Separate out the strains and aggregate the annotations per strain
;;; from original db querys

;;; (separate-strains "orig-db-genes-locus-tags-pg320.csv" "PG/PG320")
(defn separate-strains [orig-db-csv strain-dir]
  (letio [id-nm-map {"17" "locus_tag", "23" "protein_id",
                     "32" "gene", "47" "old_locus_tag"}
          tags ["gene" "old_locus_tag" "locus_tag" "protein_id"]
          cols (coll/concatv ["ACCEN" "Start" "End" "Strand"] tags)
          strainbase (fs/join panclus-base strain-dir)
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
  (let [base panclus-base
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
  (let [pg30fnas  (->> (fs/glob (fs/join panclus-base "PG/PG30/T*.fna"))
                       sort vec)
        pg320fnas (->> (fs/glob (fs/join panclus-base "PG/PG320/T*.fna"))
                       sort vec)]
    {:pg30 pg30fnas
     :pg320 pg320fnas}))

(def reffnas
  (->> (fs/glob (fs/join panclus-base "RefSeq77/N*.fna"))
       sort vec))

(def strains
  (let [nmfn #(-> % fs/basename (fs/replace-type ""))]
    {:ref77 (mapv nmfn reffnas)
     :pg30 (mapv nmfn (:pg30 pgfnas))
     :pg320 (mapv nmfn (:pg320 pgfnas))}))



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
           (let [s (if aa (ntseq2aaseq s) s)
                 [ent & x] (->> n (str/split #","))
                 lt (->> x (coll/partitionv-all 2) (into {})
                         (#(% "locus_tag")))
                 n (str ent "," lt)]
             [n (count s) (p/probs ows s)]))
         sqs))

(defn gen-strain-dists
  "Take a fasta file 'strain-fna', which contains the CDS/genes (as NT
  sqs) for one or more strains and produce a set of corresponding
  probability distributions for each CDS/gene using 'ows' kmer
  size. If 'aa' is true, first convert to ammino acid sequences.

  Returns a vector of triplets [ent len dist], one for each CDS/gene,
  where 'ent' is the entry defining the loci, 'len' is the length of
  the sq, and 'dist' is the probability distribution."
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


(defn gen-len-dists-map
  "Take a dists set (vector of triples [ent len dist] as obtained by
  gen-strain-dists) and returns a map grouping these by 'len'"
  [dists]
  (reduce (fn[S [n cnt P :as v]]
            (assoc S cnt (conj (get S cnt []) v)))
          {} dists))

(defn gen-nm-dists-map
  "Take a dists set (vector of triples [ent len dist] as obtained by
  gen-strain-dists) and returns a map grouping these by 'ent'"
  [dists]
  (reduce (fn[M [n cnt P :as v]] (assoc M n v)) {} dists))

(defn fna2center-start [strain-fna]
  (map (fn[tuple]
         {:center tuple :members #{} :name (first tuple)})
       (first (gen-strain-group-dists [strain-fna] :ows 6))))


(comment

  (def ref-dists-all (gen-strain-group-dists reffnas :ows 6))
  (def center-dists (map (fn[tuple]
                           {:center tuple :members #{} :name (first tuple)})
                         (first ref-dists-all)))
  (def ref-dists (rest ref-dists-all))
  (def refdists-map (gen-nm-dists-map ref-dists))
  (def ref-len-dists-map (gen-len-dists-map ref-dists))

  (io/with-out-writer "/store/data/PanClus/Entropy/ref-len-dists-map.clj"
    (prn ref-len-dists-map))



  (let [base "/store/data/PanClus/Entropy"
        data-nms [[center-dists "center-dists.clj"]
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


(defn ld-map-mean [ld-map]
  (->> ld-map
     (reduce (fn[[n d] [k v]]
               (let [cnt (count v)] [(+ n (* cnt k)) (+ d cnt)]))
             [0.0 0.0])
     (apply /)))

(def PGSP-index (atom 1))

(defn hybrid-center? [center-name]
  (= (str/take 4 center-name) "PGSP"))

(defn cluster-strains
  [dists center-dists
   & {:keys [len-max %-max jsdctpt] :or {len-max 5 %-max 0.99 jsdctpt 0.40}}]

  (let [dists (->> dists (apply concat))
        ref-len-dists-map (gen-len-dists-map dists)
        mean (p/mean (map second dists))
        std  (p/std-deviation (map second dists))
        mean-std (long (- mean std))
        mean+std (long (+ mean std 7))]
    (loop [center-dists center-dists
           i @PGSP-index
           clus []]
      (if (empty? center-dists)
        clus
        (let [{:keys [center members name]} (first center-dists)
              [n1 pcnt P] center
              ;;_ (prn :n1 n1 :pcnt pcnt :means mean-std mean+std)
              hybrid? (hybrid-center? name)
              ;; Adjust distribution set to use for this clu (P
              ;; 'centered') by length - Only for optimization - not
              ;; actually needed for RE
              ds (if (< mean-std pcnt mean+std)
                   (->> (range (- pcnt len-max -1)
                               (+ pcnt len-max))
                        (mapcat ref-len-dists-map))
                   (->> (range (-> pcnt (* 0.99) Math/round)
                               (-> pcnt (* 1.01) Math/round))
                        (mapcat ref-len-dists-map)))
              ;; Gen new info tuples including JSD score
              scored-tuples (sort-by
                             second
                             (vfold (fn[[n cnt Q]]
                                      [n (it/jensen-shannon P Q) cnt Q])
                                    ds))
              ;; Get new cluster memebers for center P
              new-members (coll/takev-while
                           #(-> % second (< jsdctpt))
                           scored-tuples)
              ;; Cluster for center P
              clu {:center [n1 0.0 pcnt P]
                   :members (into members new-members)
                   :name (if hybrid? name (format "PGSP_%05d" i))}]
          (recur (rest center-dists)
                 (if hybrid? i (swap! PGSP-index inc))
                 (conj clus clu)))))))

(defn cluster-leftovers
  "After a cluster pipe cycling converges there may (most likely will)
  be 'leftover points' that do not cluster with any of the clusters in
  the final set. These, as the 'dists' parameter here, will be taken
  and clustered among themselves. Centers are formed by walking
  through the 'dists' set, and the remainging points are clustered to
  each such center, and all in the resulting cluster are removed from
  'dists'. The process repeats until there are no more points. Note,
  this can stil result in singlets, which are new leftovers in the
  clutrer pipe cycling."
  [dists & {:keys [len-max %-max jsdctpt]
            :or {len-max 5 %-max 0.99 jsdctpt 0.40}}]
  (let [distmap (gen-nm-dists-map dists)
        ref-len-dists-map (gen-len-dists-map (vals distmap))
        mean (p/mean (map second dists))
        std  (p/std-deviation (map second dists))
        mean-std (long (- mean std))
        mean+std (long (+ mean std 7))]
    (loop [distmap distmap
           i @PGSP-index
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
              scored-tuples (sort-by
                             second
                             (vfold (fn[[n cnt Q]]
                                      [n (it/jensen-shannon P Q) cnt Q])
                                    ds))
              new-members (coll/takev-while
                           #(-> % second (< jsdctpt))
                           scored-tuples)
              clu {:center [n1 0.0 pcnt P]
                   :members (set new-members)
                   :name (format "PGSP_%05d" i)}
              new-distmap (apply dissoc distmap (mapv first new-members))]
          (recur new-distmap
                 (swap! PGSP-index inc)
                 (conj clus clu)))))))

(defn clean-clusters
  "Cleans the 'members' field to only have 'ent,lt' names - no counts,
  or distributions as they are not needed for cluster information."
  [clusters]
  (map (fn[clu]
         (let [members (->> clu :members (mapv first) set)]
           (assoc clu :members members)))
       clusters))

(defn make-hybrids
  "Transform the clusters in 'clusters' obtained from distribution map
  'distmap' into dist recs of form: [\"hybrid-dict-x\", cnt,
  hybrid-dict] where 'cnt' is the mean seq lens in cluster x and
  'hybrid-dict' is the minimum entropy derived hybrid dictionary /
  distribution of all distributions in cluster x"
  [clusters & {:keys [wz] :or {wz 6}}]
  (mapv (fn[clu]
         (let [{:keys [center members name]} clu
               new-members (filter coll? members)
               new+center (conj new-members center)
               old-members (filter string? members)
               dists (mapv last new+center)
               cnt (->> new+center
                        (mapv #(-> % rest second))
                        p/mean long)]
           {:center [name cnt (it/hybrid-dictionary wz dists)]
            :members (into (set old-members) (mapv first new-members))
            :name name}))
       clusters))

(defn- all-matched
  "Collect the ID set of all loci (genes/CDS) currently a member of a
  cluster. 'clusters' is a set of cluster elements, each of which has
  fields [center, members, name]. The 'members' field contains the
  loci of the cluster in either the (to be rolled into next hybrid
  center) form [n len dist] or simply n."
  [clusters]
  (->> clusters (mapcat :members)
       (map #(if (coll? %) (first %) %))))

(defn leftovers-from-clustering
  "Obtain the remaing unclustered members of 'start-dist' as those
  that are not members of any cluster in 'clusters'. 'clusters' is the
  result of cluster-strains with input 'start-dists'. Each member of
  clusters is a map with keys [center members dist] and each member of
  start-dists is standard tuple [n cnt dist]. Returns collection of
  tuples in start-dists not in any cluster."
  [clusters start-dists]
  (let [matched (all-matched clusters)]
    (vals (reduce (fn[M k]
                    (dissoc M k))
                  (gen-nm-dists-map start-dists)
                  matched))))

(defn merge-leftovers-into-clusters
  "Effectively 'leftovers' do not cluster and so are the
  start (centers) of new clusters. Clusters is the current collection
  of clusters - each a map with keys [center, members, name]. Convert
  each tuple [n cnt dist] in leftovers to such a map and add the
  result to 'clusters'. Returns new cluster coll"
  [clusters leftovers]
  (coll/concatv (mapv (fn[m]
                        (let [center (m :center)
                              [n jsd cnt dist] center]
                          (assoc m :center [n cnt dist])))
                      clusters)
                (mapv (fn[[n cnt dist]]
                        (let [i @PGSP-index]
                          (swap! PGSP-index inc)
                          {:center [n cnt dist]
                           :members #{n}
                           :name (format "PGSP_%05d" i)}))
                      leftovers)))


(defn run-clustering-pipe
  "For a single set of strain distribution records 'dists' using a
  starting set of centroid distributions (cluster centers)
  'center-dists, perform a complete cycling clustering."
  [dists center-dists
   & {:keys [len-max %-max jsdctpt diffcut]
      :or {len-max 5 %-max 0.99 jsdctpt 0.40 diffcut 5}}]
  (let [all-dists (apply concat dists)]
    (loop [clusters (cluster-strains dists center-dists :jsdctpt jsdctpt)
           leftcnt (count all-dists)]
      #_(prn "         (first clusters)" (first clusters))
      (let [hybrids (make-hybrids clusters)
            leftovers (leftovers-from-clustering clusters all-dists)
            cur-leftcnt (count leftovers)
            nextclus (cluster-strains [leftovers] hybrids)]
        (println :diff (- leftcnt cur-leftcnt))
        (if (< (- leftcnt cur-leftcnt) diffcut)
          nextclus
          (recur nextclus, cur-leftcnt))))))

;;; (map (fn[tuple]
;;;        {:center tuple :members #{} :name (first tuple)})
;;;      (gen-strain-group-dists [start-center-fna] :ows ows))
(defn run-strain-clustering
  ""
  [strain-fnas center-start
   & {:keys [len-max %-max jsdctpt diffcut chunk-size]
      :or {len-max 5 %-max 0.99 jsdctpt 0.40 diffcut 5 chunk-size 10}}]
  (let [alpha :aa
        ows (ows alpha)]
    (loop [strain-fnas strain-fnas
           clusters center-start]
      (if (empty? strain-fnas)
        clusters
        (let [chunk-fnas (take chunk-size strain-fnas)
              chunk-dists (gen-strain-group-dists chunk-fnas :ows ows)
              chunk-clusters (run-clustering-pipe
                              chunk-dists clusters
                              :diffcut diffcut :jsdctpt jsdctpt)
              leftover-dists (leftovers-from-clustering
                              chunk-clusters (apply concat chunk-dists))
              leftover-clus (cluster-leftovers
                             leftover-dists :jsdctpt jsdctpt)
              leftover-hybs (make-hybrids
                             (filter #(> (count (% :members)) 1) leftover-clus))
              _ (println "After leftover-hybs")
              clusters (merge-leftovers-into-clusters
                        chunk-clusters leftover-dists)]
          (println "At recur point")
          (recur (drop chunk-size strain-fnas)
                 clusters))))))



(comment

(def ref77-SP-clusters
  (time (run-strain-clustering
         (rest reffnas) (fna2center-start (first reffnas))
         :jsdctpt 0.6 :chunk-size 20)))

(def ref77+PG30-SP-clusters
  (time (run-strain-clustering
         (pgfnas :pg30) ref77-SP-clusters
         :jsdctpt 0.6 :chunk-size 30)))

(def ref77+PG30-SP-clusters
  (time (run-strain-clustering
         (concat (rest reffnas) (pgfnas :pg30))
         (fna2center-start (first reffnas))
         :jsdctpt 0.6 :chunk-size 50)))
(def ref77+PG30-SP-clusters
  (->> ref77+PG30-SP-clusters (filter #(> (count (% :members)) 1))))


(def leftovers-ref77pg30
  (->> ref77+PG30-SP-clusters (filter #(= (count (% :members)) 1))
       (map (fn[m] (let [[nm jsd cnt dist] (m :center)
                        n (first (m :members))
                        [cnt dist] (if dist [cnt dist] [jsd cnt])]
                    [n cnt dist])))))

(def leftover-clus-ref77+pg30
  (cluster-leftovers leftovers-ref77pg30 :jsdctpt 0.6))
(def jsa (let [x (make-hybrids
                  (->> leftover-clus-ref77+pg30
                       (filter #(> (count (% :members)) 1))))]
           (run-clustering-pipe [leftovers-ref77pg30] x)))



(->> ref77-SP-clusters (map :members) (filter #(> (count %) 19)) count)
1179
(->> ref77-SP-clusters (map :members) (filter #(> (count %) 18)) count)
1535

(->> ref77+PG30-SP-clusters (map :members) (filter #(> (count %) 47)) count)
1376
(->> ref77+PG30-SP-clusters (map :members) (filter #(> (count %) 48)) count)
1297
(->> ref77+PG30-SP-clusters (map :members) (filter #(> (count %) 49)) count)
1029



(def ref-10-clusters-tst
  (time (cluster-strains (take 10 ref-dists) center-dists :jsdctpt 0.40)))
(def all-matched-tst
  (->> ref-10-clusters-tst (mapcat :members) (map first)))
(def hybrid-tst (make-hybrids ref-10-clusters-tst))
(def leftover-clus
  (cluster-strains
   [(leftovers-from-clustering
     ref-10-clusters-tst (->> ref-dists (take 10) (apply concat)))]
   hybrid-tst))
(io/with-out-writer
  "/store/data/PanClus/Entropy/ref-clusters-040.clj"
  (prn ref-clusters-040))

(count (leftovers-from-clustering
        ref-10-clusters-tst (->> ref-dists (take 10) (apply concat))))
(count (leftovers-from-clustering
        leftover-clus (->> ref-dists (take 10) (apply concat))))
(count (leftovers-from-clustering
        leftover-clus-2 (->> ref-dists (take 10) (apply concat))))


(->> ref-10-clusters-tst (map :members) (map count) m/sum)
18895
(->> leftover-clus (map :members) (map count) m/sum)
21637
(->> leftover-clus (map :members) (map count) (filter #(> % 0)) count)
2216
(->> ref-10-clusters-tst (map :members) (map count) (filter #(> % 0)) count)
2043



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


(cluster-strains ref-dists (take 5 (coll/dropv 1000 center-dists)))
(cluster-strains ref-dists (take 5 (coll/dropv 2000 center-dists)))


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
(->> center-dists (take-while #(< (second %) 60)) count)
120
(->> center-dists (take-while #(< (second %) 504)) count)

(def jsadists
  (conj (map refdists-map (butlast (map first (first *1))))
        (first (coll/dropv 1000 center-dists))))

(p/mean (map second (map refdists-map (map first (rest jsadists)))))

(def hybdist
  ["hybrid-dist-test" 200 (it/hybrid-dictionary 0 (map last jsadists))])

(count (last hybdist))
236
(map count (map last (map refdists-map (map first (rest jsadists)))))
(195 195 195 195 195 195 195 195 195 195 195 195 194)

(cluster-strains ref-dists [hybdist])
)
