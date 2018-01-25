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

(defn member->entry
  "Return entry part of member id 'member' which is entry,locus-tag"
  [member]
  (->> member (str/split #",") first))

(defn member->locus
  "Return locus-tag part of member id 'member' which is entry,locus-tag"
  [member]
  (->> member (str/split #",") second))

(defn member->strain
  "Return strain part of entry of member id 'member' which is entry,locus-tag"
  [member]
  (->> member (str/split #",") first bufiles/entry-parts first))


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


(defn running-mean [{:keys [sm d m] :as M} n]
  (let [sm (+ sm n), d (inc d)]
    (assoc M :sm sm :d d :m (-> sm (/ d) double))
    #_{:sm sm, :d d, :m (-> sm (/ d) double)}))


(def PGSP-index (atom 0))

(defn init-clusters
  "Initialize a map of clusters based on content of vector/map
  coll. If coll is a map, assumed to already be a correctly formatted
  map of clusters and simply returns this.

  If (vector? coll), each element must have form [n cnt P], where n is
  the 'entry,locus_tag' name, cnt is the length of sequence, and P the
  kmer distribution.

  Returns a map with entries [k m], with k a PGSP locus tag formed
  using 'pgsp-start' integer as beginning number, and m a map of the
  form:

  {:center [k cnt Q], :dists [Q] :jsds [n 0.0]
   :cnt {:sm cnt, :d 1, :m cnt}
   :members #{n} :name k}"

  [coll pgsp-start]
  (if (map? coll)
    coll
    (do
      (swap! PGSP-index (constantly pgsp-start))
      (reduce (fn[M [nm cnt Q]]
                (let [pgnm (format "PGSP_%05d" (swap! PGSP-index inc))]
                  (assoc M pgnm {:center [pgnm cnt Q]
                                 :dists [Q]
                                 :jsds [[nm 0.0]]
                                 :cnt {:sm cnt, :d 1, :m cnt}
                                 :members #{nm}
                                 :name pgnm})))
              {} coll))))


(defn cluster-strains
  ""
  [dists start-clusters
   & {:keys [len-max %-max jsdctpt pgsp-start]
      :or {len-max 5 %-max 0.99 jsdctpt 0.40 pgsp-start @PGSP-index}}]

  (let [dists (->> dists (apply concat))
        clus (init-clusters start-clusters pgsp-start)
        center-dists (mapv #(% :center) (vals clus))
        len-distmap (gen-len-dists-map center-dists)]
    (loop [dists dists
           clus clus]
      (if (empty? dists)
        clus
        (let [[n1 pcnt P] (first dists)
              ;;_ (prn :n1 n1 :pcnt pcnt :means mean-std mean+std)

              ;; Adjust clus distribution set to use, for current
              ;; dists point, by length - Only for optimization - not
              ;; actually needed for RE
              ds (->> (range (-> pcnt (* 0.80) Math/round)
                             (-> pcnt (* 1.20) Math/round))
                      (mapcat len-distmap))
              ;; Gen scores of current dists point to centers in 'ds'
              ;; set. Pick center with lowest score.
              best-center (->> ds
                               (vfold (fn[[n cnt Q]]
                                        [n (it/jensen-shannon P Q) cnt Q]))
                               (sort-by second)
                               first)
              ;; If best-center score <= jsdctpt, add current point to
              ;; the cluster for this center
              clu (when (and best-center (<= (second best-center) jsdctpt))
                    (let [[nm jsd cnt Q] best-center
                          curclus (get clus nm)
                          members (:members curclus)
                          new-dists (:dists curclus)
                          jsds (conj (curclus :jsds) [n1 jsd])
                          new-cnt (-> curclus :cnt (running-mean pcnt))]
                      (-> curclus
                          (assoc-in [:dists] (conj new-dists P))
                          (assoc-in  [:members] (conj members n1))
                          (assoc-in [:jsds] jsds)
                          (assoc-in [:cnt] new-cnt))))
              new-clus (if clu (assoc clus (clu :name) clu) clus)]
          #_(println :new-distmap (count new-distmap))
          (recur (rest dists)
                 new-clus))))))


(defn make-hybrids
  "Transform the clusters in 'clusters' as obtained from
  'cluster-strains' to have updated information (hybrid) centroid
  distributions suitable for new rounds of clustering."
  [clusters & {:keys [wz] :or {wz 6}}]
  (reduce (fn[M clu]
            (let [{:keys [center members dists jsds cnt name]} clu
                  Q (it/hybrid-dictionary wz dists)]
              (assoc M name {:center [name (-> cnt :m long) Q]
                             :dists [Q]
                             :jsds jsds
                             :cnt cnt
                             :members members
                             :name name})))
          {} (if (map? clusters) (vals clusters) clusters)))


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
  (let [clusters (if (map? clusters) (vals clusters) clusters)
        matched (all-matched clusters)]
    (vals (reduce (fn[M k]
                    (dissoc M k))
                  (gen-nm-dists-map start-dists)
                  matched))))


(defn run-clustering-cycle
  "For a single set of strain distribution records 'dists' using a
  starting set of centroid distributions (cluster centers)
  'center-dists, perform a complete cycling clustering."
  [dists center-dists
   & {:keys [len-max %-max jsdctpt diffcut pgsp-start]
      :or {len-max 5 %-max 0.99 jsdctpt 0.40 diffcut 0}}]
  (let [all-dists (apply concat dists)]
    (loop [clusters (cluster-strains
                     dists center-dists
                     :len-max len-max :%-max %-max
                     :jsdctpt jsdctpt :pgsp-start pgsp-start)
           leftcnt (count all-dists)]
      #_(prn "         (first clusters)" (first clusters))
      (let [hybrids (make-hybrids clusters)
            leftovers (leftovers-from-clustering clusters all-dists)
            cur-leftcnt (count leftovers)
            nextclus (cluster-strains [leftovers] hybrids)]
        (println :diff (- leftcnt cur-leftcnt))
        (if (<= (- leftcnt cur-leftcnt) diffcut)
          nextclus
          (recur nextclus, cur-leftcnt))))))


(defn merge-clusters
  "Effectively 'leftovers' do not cluster and so are the
  start (centers) of new clusters. Clusters is the current collection
  of clusters - each a map with keys [center, members, name]. Convert
  each tuple [n cnt dist] in leftovers to such a map and add the
  result to 'clusters'. Returns new cluster coll"
  ([clus1 clus2]
   (merge clus1 clus2))
  ([clus1 clus2 & clusters]
   (apply merge clus1 clus2 clusters)))


(defn get-strain-id [dist-id]
  (->> dist-id (str/split #",") first (str/split #"/") first))

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
  (if (empty? dists)
    {}
    (let [strain-map (group-by (fn[v] (->> v first get-strain-id)) dists)
          strain (->> strain-map
                      (mapv (fn[[k v]] [k (count v)]))
                      (sort-by second) ffirst)
          new-starts (strain-map strain)
          nxt-dists (-> strain-map (dissoc strain) vals)
          new-clus  (run-clustering-cycle
                     nxt-dists new-starts
                     :len-max len-max :%-max %-max
                     :jsdctpt jsdctpt :pgsp-start @PGSP-index)]
      (merge-clusters
       new-clus
       (cluster-leftovers (leftovers-from-clustering new-clus dists)
                          :len-max len-max :%-max %-max
                          :jsdctpt jsdctpt)))))


(defn clustering->singlet-dists
  [clustering]
  (->> clustering vals
       (filter #(-> % :members count (= 1)))
       (mapv (fn[m]
               (let [[nm cnt dist] (m :center)
                     name   (-> m :members first)]
                 [name cnt dist])))))

(defn clustering-minus-singlets
  [clustering]
  (->> clustering
       (filter (fn[[n m]]
                 (-> m :members count (> 1))))
       (into {})))

(defn final-pass
  [clustering]
  (let [singlet-dists (clustering->singlet-dists clustering)
        start-centers (clustering-minus-singlets clustering)
        pgsp-num (->> start-centers keys sort last
                      (str/split #"_") second Integer.)
        newclus (run-clustering-cycle
                 [singlet-dists] start-centers
                 :pgsp-start pgsp-num)]
    (merge-clusters
     newclus
     (cluster-leftovers (leftovers-from-clustering newclus singlet-dists)))))


(defn run-strain-clustering
  ""
  [strain-fnas centers-start
   & {:keys [len-max %-max jsdctpt diffcut chunk-size pgsp-start]
      :or {len-max 5 %-max 0.99 jsdctpt 0.40 diffcut 0 chunk-size 10}}]
  (println :let)
  (let [alpha :aa
        ows (ows alpha)]
    (println :at-loop)
    (loop [pgsp-num pgsp-start
           strain-fnas strain-fnas
           clusters centers-start]
      (if (empty? strain-fnas)
        (final-pass clusters)
        (let [chunk-fnas (take chunk-size strain-fnas)
              _ (println :fnas)
              chunk-dists (gen-strain-group-dists chunk-fnas :ows ows)
              _ (println :dists)

              chunk-clusters (run-clustering-cycle
                              chunk-dists clusters
                              :len-max len-max :%-max %-max
                              :jsdctpt jsdctpt :pgsp-start pgsp-num
                              :diffcut diffcut)

              leftover-clus (cluster-leftovers
                             (leftovers-from-clustering
                              chunk-clusters (apply concat chunk-dists))
                             :len-max len-max :%-max %-max :jsdctpt jsdctpt)]
          (println "At recur point")
          (recur @PGSP-index
                 (drop chunk-size strain-fnas)
                 (merge-clusters chunk-clusters leftover-clus)))))))


(defn clean-clusters
  "Cleans the 'members' field to only have 'ent,lt' names - no counts,
  or distributions as they are not needed for cluster information."
  [clusters]
  (map (fn[clu]
         (let [members (->> clu :members (mapv first) set)]
           (assoc clu :members members)))
       clusters))


(defn clus-jsd-data
  [clusters & {:keys [cnt] :or {cnt 1}}]
  (->> clusters
       (filter (fn[[k m]] (-> m :members count (>  cnt))))
       vals
       (mapv (fn[m]
               (let [jsds (->> m :jsds (sort-by second) (mapv second))
                     mjsd (p/mean jsds)
                     mxjsd (last jsds)
                     vjsd (p/variance jsds)
                     stdjsd (p/std-deviation jsds)
                     n1sdev (- mjsd stdjsd)
                     p1sdev (+ mjsd stdjsd)
                     sdev-jsds (filter #(<= n1sdev % p1sdev) jsds)
                     less+1sdev-jsds (filter #(<= % p1sdev) jsds)]
                 {:mean mjsd :mx mxjsd :var vjsd :sdev stdjsd
                  :-1sdev n1sdev :+1sdev p1sdev
                  :jsds jsds
                  :within-1sdev sdev-jsds
                  :less+1sdev less+1sdev-jsds})))))

(defn clus-counts
  [clusters & {:keys [cnt sdev-grp cut%]
               :or {cnt 1, sdev-grp :within-1sdev cut% 0.8}}]
  (let [jsdata (clus-jsd-data clusters :cnt cnt)
        total (->> jsdata
                   (filter #(>= (-> % sdev-grp count
                                    (/ (-> % :jsds count))
                                    double)
                                cut%))
                   count)
        percent (-> total (/ (count jsdata)) double roundit)]
    [total percent]))

(defn clus-jsd-dists
  [clusters]
  )



(comment

(let [x (leftovers-from-clustering
         newclu
         (->> ref-dists (take 2) (apply concat)))
      strain-map (group-by (fn[v] (->> v first get-strain-id)) x)
      strain (->> strain-map
                  (mapv (fn[[k v]] [k (count v)]))
                  (sort-by second) ffirst)
      new-starts (strain-map strain)
      nxt-dists (-> strain-map (dissoc strain) vals)]
  [(count new-starts) (count nxt-dists) (->> nxt-dists (apply concat) count)])

(def jsa
  (let [x (leftovers-from-clustering
           newclu
           (->> ref-dists (take 2) (apply concat)))
        left-clus (cluster-leftovers x)]
    left-clus))

(def newclu
  (cluster-strains (take 2 ref-dists) (first ref-dists-all) :pgsp-start 0))




(def ref77-SP-clusters
  (time (run-strain-clustering
         (rest reffnas)
         (-> reffnas first (gen-strain-dists :ows (ows :aa)))
         :pgsp-start 0 :jsdctpt 0.4 :chunk-size 2 :diffcut 0)))
"Elapsed time: 8545367.909075 msecs"
(def ref77-SP-clusters
  (-> (pams/get-params :panclus-base)
      (fs/join "Entropy/ref77-SP-clusters.clj")
      slurp read-string))

(def ref77pg30-fut
  (future (run-strain-clustering
           (pgfnas :pg30)
           ref77-SP-clusters
           :pgsp-start @PGSP-index :jsdctpt 0.4 :chunk-size 2 :diffcut 0)))

(def ref77+PG30-SP-clusters
  (time (run-strain-clustering
         (coll/concatv (rest reffnas) (pgfnas :pg30))
         #_(-> reffnas first (gen-strain-dists :ows (ows :aa)))
         ref77-SP-clusters
         :pgsp-start 0 :jsdctpt 0.4 :chunk-size 2 :diffcut 0)))

(def ref77pg30-nosinglets (clustering-minus-singlets ref77+PG30-SP-clusters))
(def singlets (clustering->singlet-dists ref77+PG30-SP-clusters))


(def ref77+PG30-SP-clusters
  (-> (pams/get-params :panclus-base)
      (fs/join "Entropy/ref77+PG30-SP-clusters.clj")
      slurp read-string))
(def ref77pg350-fut
  (future (run-strain-clustering
           (pgfnas :pg320)
           ref77+PG30-SP-clusters
           :pgsp-start @PGSP-index :jsdctpt 0.55 :chunk-size 2 :diffcut 0)))







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
