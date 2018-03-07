(ns panclus.stats
  [:require
   [clojure.data.csv :as csv]
   [clojure.string :as cljstr]
   [clojure.set :as set]

   [aerial.fs :as fs]

   [aerial.utils.string :as str]
   [aerial.utils.io :refer [letio] :as io]
   [aerial.utils.coll :refer [vfold] :as coll]
   [aerial.utils.math :as m]
   [aerial.utils.math.probs-stats :as p]
   [aerial.utils.math.infoth :as it]
   [aerial.utils.math.combinatorics :as cmb]

   [aerial.bio.utils.files :as bufiles]

   [panclus.params :as pams]
   [panclus.clustering
    :refer :all
            #_[member->entry
            member->strain
            member->locus
            member->sq
            member->dist
            member->clus-input
            ntseq2aaseq
            id-line->member
            roundit
            PGSP-index
            make-hybrids
            get-dist-cnt
            get-dist
            running-mean
            run-strain-clustering
            pair-jsd
            next-set
            tran-reduce
            merged-set->cluster
            merged-sets->clusters
            #_ref77-SP-clusters
            #_ref77+PG30-SP-clusters
            ref77pg350-fut
            ]]]
  )



(defn clus-jsd-data
  [clusters & {:keys [cnt] :or {cnt 1}}]
  (->> clusters
       (filter (fn[[k m]] (-> m :members count (>  cnt))))
       vals
       (mapv (fn[m]
               (let [jsds (->> m :jsds (sort-by second) (mapv second))
                     outjsds (->> m :outjsds)
                     memcnt (->> m :members count)
                     pgnm (m :name)
                     mjsd (p/mean jsds)
                     mxjsd (last jsds)
                     vjsd (p/variance jsds)
                     stdjsd (p/std-deviation jsds)
                     n1sdev (- mjsd stdjsd)
                     p1sdev (+ mjsd stdjsd)
                     sdev-jsds (filter #(<= n1sdev % p1sdev) jsds)
                     less+1sdev-jsds (filter #(<= % p1sdev) jsds)]
                 {:pgnm pgnm :memcnt memcnt
                  :mean mjsd :mx mxjsd :var vjsd :sdev stdjsd
                  :-1sdev n1sdev :+1sdev p1sdev
                  :jsds jsds
                  :outjsds outjsds
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

(clus-counts ref77-SP-clusters :sdev-grp :less+1sdev :cut% 0.8)
(clus-counts ref77-SP-clusters :sdev-grp :within-1sdev :cut% 0.8)

(clus-counts ref77+PG30-SP-clusters :sdev-grp :less+1sdev :cut% 0.8)
(clus-counts ref77+PG30-SP-clusters :sdev-grp :within-1sdev :cut% 0.8)

(clus-counts @ref77pg350-fut :sdev-grp :less+1sdev :cut% 0.8)
(clus-counts @ref77pg350-fut :sdev-grp :within-1sdev :cut% 0.8)

(clus-counts (@DBG :current-merged-clustering)
             :sdev-grp :less+1sdev :cut% 0.8) ; 89%
(clus-counts (@DBG :current-merged-clustering) ; 65%
             :sdev-grp :within-1sdev :cut% 0.8)




;;; cluster seq and kmcnt distributions
(io/with-out-writer "/store/data/PanClus/Stats/ref77-len-dist.clj"
  (prn (->>  ref77-SP-clusters vals (mapv #(-> % :center second))
             (group-by identity)
             (mapv (fn[[len elts]] {:len len :cnt (count elts)})))))

(io/with-out-writer "/store/data/PanClus/Stats/ref77-kmcnt-dist.clj"
  (prn (->>  ref77-SP-clusters vals (mapv #(-> % :center last))
             (group-by #(count %))
             (mapv (fn[[kmcnt elts]] {:kmcnt kmcnt :cnt (count elts)})))))


(io/with-out-writer "/store/data/PanClus/Stats/ref77-pg30-len-dist.clj"
  (prn (->>  ref77+PG30-SP-clusters vals (mapv #(-> % :center second))
             (group-by identity)
             (mapv (fn[[len elts]] {:len len :cnt (count elts)})))))

(io/with-out-writer "/store/data/PanClus/Stats/ref77-pg30-kmcnt-dist.clj"
  (prn (->>  ref77+PG30-SP-clusters vals (mapv #(-> % :center last))
             (group-by #(count %))
             (mapv (fn[[kmcnt elts]] {:kmcnt kmcnt :cnt (count elts)})))))


(io/with-out-writer "/store/data/PanClus/Stats/ref77-pg350-len-dist.clj"
  (prn (->>  @ref77pg350-fut vals (mapv #(-> % :center second))
             (group-by identity)
             (mapv (fn[[len elts]] {:len len :cnt (count elts)})))))

(io/with-out-writer "/store/data/PanClus/Stats/ref77-pg350-kmcnt-dist.clj"
  (prn (->>  @ref77pg350-fut vals (mapv #(-> % :center last))
             (group-by #(count %))
             (mapv (fn[[kmcnt elts]] {:kmcnt kmcnt :cnt (count elts)})))))




;;; Cluster distributions
(io/with-out-writer "/store/data/PanClus/Stats/ref77-clustering-dist.clj"
  (prn (->> ref77-SP-clusters vals
            (group-by #(count (% :members)))
            (map (fn[[k m]] [k (count m)]))
            (sort-by first)
            (mapv (fn[[sz ct]] {:sz sz :cnt ct})))))


(io/with-out-writer "/store/data/PanClus/Stats/ref77-pg30-clustering-dist.clj"
  (prn (->> ref77+PG30-SP-clusters vals
            (group-by #(count (% :members)))
            (map (fn[[k m]] [k (count m)]))
            (sort-by first)
            (mapv (fn[[sz ct]] {:sz sz :cnt ct})))))


(io/with-out-writer "/store/data/PanClus/Stats/ref77-pg350-clustering-dist.clj"
  (prn (->> @ref77pg350-fut vals
            (group-by #(count (% :members)))
            (map (fn[[k m]] [k (count m)]))
            (sort-by first)
            (mapv (fn[[sz ct]] {:sz sz :cnt ct})))))

(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg350-merged-clustering-dist.clj"
  (prn (->> (@DBG :current-merged-clustering) vals
            (group-by #(count (% :members)))
            (map (fn[[k m]] [k (count m)]))
            (sort-by first)
            (mapv (fn[[sz ct]] {:sz sz :cnt ct})))))





;;; Cluster radii - average JSD to center
(io/with-out-writer "/store/data/PanClus/Stats/ref77-clus-radii.clj"
  (prn (->> ref77-SP-clusters
            clus-jsd-data
            (group-by #(-> % :mean (roundit :places 2)))
            (map (fn[[k m]] [k (count m)]))
            (sort-by first)
            (mapv (fn[[jsd ct]] {:jsd jsd :cnt ct})))))


(io/with-out-writer "/store/data/PanClus/Stats/ref77-pg30-clus-radii.clj"
  (prn (->> ref77+PG30-SP-clusters
            clus-jsd-data
            (group-by #(-> % :mean (roundit :places 2)))
            (map (fn[[k m]] [k (count m)]))
            (sort-by first)
            (mapv (fn[[jsd ct]] {:jsd jsd :cnt ct})))))


(io/with-out-writer "/store/data/PanClus/Stats/ref77-pg350-clus-radii.clj"
  (prn (->> @ref77pg350-fut
            clus-jsd-data
            (group-by #(-> % :mean (roundit :places 2)))
            (map (fn[[k m]] [k (count m)]))
            (sort-by first)
            (mapv (fn[[jsd ct]] {:jsd jsd :cnt ct})))))

(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg350-merged-clus-radii.clj"
  (prn (->> ref77-pg350-with-merged-clusters
            clus-jsd-data
            (group-by #(-> % :mean (roundit :places 2)))
            (map (fn[[k m]] [k (count m)]))
            (sort-by first)
            (mapv (fn[[jsd ct]] {:jsd jsd :cnt ct})))))




;;; Cluster max pt - JSD to center
(io/with-out-writer "/store/data/PanClus/Stats/ref77-max-outliers.clj"
  (prn (->> ref77-SP-clusters
            (#(clus-jsd-data % :cnt 1))
            (group-by #(-> % :mx (roundit :places 2)))
            (map (fn[[k m]] [k (count m)]))
            (sort-by first)
            (mapv (fn[[mx ct]] {:mx mx :cnt ct}))
            #_(filter #(<= (% :mx) 0.3))
            #_(mapv :cnt)
            #_m/sum))) ==> 92% <= .3

(io/with-out-writer "/store/data/PanClus/Stats/ref77-pg30-max-outliers.clj"
  (prn (->> ref77+PG30-SP-clusters
            (#(clus-jsd-data % :cnt 1))
            (group-by #(-> % :mx (roundit :places 2)))
            (map (fn[[k m]] [k (count m)]))
            (sort-by first)
            (mapv (fn[[mx ct]] {:mx mx :cnt ct}))
            #_(filter #(<= (% :mx) 0.3))
            #_(mapv :cnt)
            #_m/sum))) ==> 90% <= .3

(io/with-out-writer "/store/data/PanClus/Stats/ref77-pg350-max-outliers.clj"
  (prn (->> @ref77pg350-fut
            clus-jsd-data
            (group-by #(-> % :mx (roundit :places 2)))
            (map (fn[[k m]] [k (count m)]))
            (sort-by first)
            (mapv (fn[[mx ct]] {:mx mx :cnt ct}))
            #_(filter #(<= (% :mx) 0.3))
            #_(mapv :cnt)
            #_m/sum))) ==> 86% <= .3

(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg350-merged-max-outliers.clj"
  (prn (->> ref77-pg350-with-merged-clusters
            clus-jsd-data
            (group-by #(-> % :mx (roundit :places 2)))
            (map (fn[[k m]] [k (count m)]))
            (sort-by first)
            (mapv (fn[[mx ct]] {:mx mx :cnt ct}))
            #_(filter #(<= (% :mx) 0.3))
            #_(mapv :cnt)
            #_m/sum))) ==> 86% <= .3

;; ==> 86% <= .3 / 99% <= .4 (54 > .4)





;;; Cluster JSD standard deviations

(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-jsd-sdevs.clj"
  (prn (->> (clus-jsd-data ref77-SP-clusters :cnt 1)
            (mapv #(vector (-> % :mean roundit) (-> % :sdev roundit)))
            (sort-by second)
            (group-by #(-> % second))
            (mapv (fn[[sdev v]] {:sdev sdev :cnt (count v)})))))

(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg30-jsd-sdevs.clj"
  (prn (->> (clus-jsd-data ref77+PG30-SP-clusters :cnt 1)
            (mapv #(vector (-> % :mean roundit) (-> % :sdev roundit)))
            (sort-by second)
            (group-by #(-> % second))
            (mapv (fn[[sdev v]] {:sdev sdev :cnt (count v)})))))

(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg350-jsd-sdevs.clj"
  (prn (->> (clus-jsd-data @ref77pg350-fut :cnt 1)
            (mapv #(vector (-> % :mean roundit) (-> % :sdev roundit)))
            (sort-by second)
            (group-by #(-> % second))
            (mapv (fn[[sdev v]] {:sdev sdev :cnt (count v)})))))

(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg350-merged-jsd-sdevs.clj"
  (prn (->> (clus-jsd-data ref77-pg350-with-merged-clusters :cnt 1)
            (mapv #(vector (-> % :mean roundit) (-> % :sdev roundit)))
            (sort-by second)
            (group-by #(-> % second))
            (mapv (fn[[sdev v]] {:sdev sdev :cnt (count v)})))))



(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg350-10-371-jsd-sdevs.clj"
  (prn (->> (clus-jsd-data @ref77pg350-fut :cnt 9)
            (filter #(< (-> % :jsds count) 372))
            (mapv #(vector (-> % :mean roundit) (-> % :sdev roundit)))
            (sort-by second)
            (group-by #(-> % second))
            (mapv (fn[[sdev v]] {:sdev sdev :cnt (count v)})))))

(io/with-out-writer
  "/store/data/PanClus/Stats/cluster-jsd-sdevs-371-only.clj"
  (prn (->> (clus-jsd-data @ref77pg350-fut :cnt 370)
            (filter #(< (-> % :jsds count) 372))
            (mapv #(vector (-> % :mean roundit) (-> % :sdev roundit)))
            (sort-by second)
            #_(coll/drop-until #(-> % second (> 0.1)))
            (group-by #(-> % second))
            (mapv (fn[[sdev v]] {:sdev sdev :cnt (count v)})))))
(io/with-out-writer
  "/store/data/PanClus/Stats/last-merged-jsd-sdevs-371-only.clj"
  (prn (->> (clus-jsd-data ref77-pg350-with-merged-clusters :cnt 370)
            (filter #(< (-> % :jsds count) 372))
            (mapv #(vector (-> % :mean roundit) (-> % :sdev roundit)))
            (sort-by second)
            #_(coll/drop-until #(-> % second (> 0.1)))
            (group-by #(-> % second))
            (mapv (fn[[sdev v]] {:sdev sdev :cnt (count v)})))))






;;; Percent within % of 1 stdev
(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-percent-clusters-1stdev.clj"
  (prn (->> (range 0.1 0.9 0.1) (map roundit)
            (mapv #(vector % (clus-counts ref77-SP-clusters
                                          :sdev-grp :within-1sdev :cut% %)))
            (mapv (fn[[per [cnt cper]]] {:cut per :cnt cnt :per cper})))))

(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg30-percent-clusters-1stdev.clj"
  (prn (->> (range 0.1 0.9 0.1) (map roundit)
            (mapv #(vector % (clus-counts ref77+PG30-SP-clusters
                                          :sdev-grp :within-1sdev :cut% %)))
            (mapv (fn[[per [cnt cper]]] {:cut per :cnt cnt :per cper})))))

(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg350-percent-clusters-1stdev.clj"
  (prn (->> (range 0.1 0.9 0.1) (map roundit)
            (mapv #(vector % (clus-counts @ref77pg350-fut
                                          :sdev-grp :within-1sdev :cut% %)))
            (mapv (fn[[per [cnt cper]]] {:cut per :cnt cnt :per cper})))))

(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg350-merged-percent-clusters-1stdev.clj"
  (prn (->> (range 0.1 0.9 0.1) (map roundit)
            (mapv #(vector % (clus-counts ref77-pg350-with-merged-clusters
                                          :sdev-grp :within-1sdev :cut% %)))
            (mapv (fn[[per [cnt cper]]] {:cut per :cnt cnt :per cper})))))






;;; Percent less than +1 stdev
(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-percent-clusters-less+1stdev.clj"
  (prn (->> (range 0.1 0.9 0.1) (map roundit)
            (mapv #(vector % (clus-counts ref77-SP-clusters
                                          :sdev-grp :less+1sdev :cut% %)))
            (mapv (fn[[per [cnt cper]]] {:cut per :cnt cnt :per cper})))))

(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg30-percent-clusters-less+1stdev.clj"
  (prn (->> (range 0.1 0.9 0.1) (map roundit)
            (mapv #(vector % (clus-counts ref77+PG30-SP-clusters
                                          :sdev-grp :less+1sdev :cut% %)))
            (mapv (fn[[per [cnt cper]]] {:cut per :cnt cnt :per cper})))))

(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg350-percent-clusters-less+1stdev.clj"
  (prn (->> (range 0.1 0.9 0.1) (map roundit)
            (mapv #(vector % (clus-counts @ref77pg350-fut
                                          :sdev-grp :less+1sdev :cut% %)))
            (mapv (fn[[per [cnt cper]]] {:cut per :cnt cnt :per cper})))))

(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg350-merged-percent-clusters-less+1stdev.clj"
  (prn (->> (range 0.1 0.9 0.1) (map roundit)
            (mapv #(vector % (clus-counts ref77-pg350-with-merged-clusters
                                          :sdev-grp :less+1sdev :cut% %)))
            (mapv (fn[[per [cnt cper]]] {:cut per :cnt cnt :per cper})))))






(->> ref77-pg350-with-merged-clusters vals
     (filter #(-> % :members count (= 1)))
     (map #(-> % :members first))
     (group-by member->strain)
     vals (mapv count))

(->> ref77-pg350-with-merged-clusters vals (filter #(-> % :members count (= 1)))
     (map #(-> % :members first))
     (group-by member->strain)
     (mapv (fn[[lt genes]] [lt (count genes)])))

(->> ref77-pg350-with-merged-clusters vals (filter #(-> % :members count (= 1)))
     (map #(-> % :members first))
     (group-by member->strain)
     (mapv (fn[[lt genes]] [lt (count genes)]))
     #_(sort-by second)
     (mapv second)
     #_p/mean
     p/median)




;;; Center distance data and distributions

(defn singlet-clu->distpt [clu]
  (let [nm (-> clu :members first)
        [_ cnt dist] (clu :center)]
    [nm cnt dist]))

(defn next-point-dist [clusters mode]
  (if (= mode :center-center)
    (-> clusters first last)
    (assert :NYI "only :center-center supported at this point")))

(defn compare-clusters [clusters & {:keys [mode] :or {mode :center-center}}]
  (loop [comp-dist {}
         clusters (mapv #(vector (% :name) (-> % :center last)) clusters)]
    (if (empty? (rest clusters))
      comp-dist
      (let [P (next-point-dist clusters mode)
            Pnm (-> clusters first first)
            dists (rest clusters)
            [Qnm score] (->> dists
                             (vfold (fn[[nm Q]]
                                      [nm (roundit (it/jensen-shannon P Q))]))
                             (sort-by second) first)]
        (println (format "(count clusters): %s" (count clusters)))
        (recur (assoc comp-dist score
                      (conj (get comp-dist score []) [Pnm Qnm]))
               (rest clusters))))))



;;; Singlets only
(def center-jsd-dist
  (->> ref77-SP-clusters
       (filter (fn[[k m]] (-> m :members count (= 1))))
       vals
       compare-clusters))

(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-singlet-center-distance-dist.clj"
  (prn (->> center-jsd-dist (sort-by first)
            (mapv (fn[[jsd members]] {:jsd jsd :cnt (count members)})))))
(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg30-singlet-center-distance-dist.clj"
  (prn (->> center-jsd-dist (sort-by first)
            (mapv (fn[[jsd members]] {:jsd jsd :cnt (count members)})))))


(def all-singlet-center-dist center-jsd-dist)
(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg350-singlet-center-distance-data.clj"
  (prn all-singlet-center-dist))
(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg350-singlet-center-distance-dist.clj"
  (prn (->> center-jsd-dist (sort-by first)
            (mapv (fn[[jsd members]] {:jsd jsd :cnt (count members)})))))



;;; Full clustering - singlets and multis
(def center-jsd-dist
  (->> ref77-SP-clusters
       #_(filter (fn[[k m]] (-> m :members count (> 1))))
       vals
       compare-clusters))
(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-center-distance-data.clj"
  (prn center-jsd-dist))

(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-center-distance-dist.clj"
  (prn (->> center-jsd-dist (sort-by first)
            (mapv (fn[[jsd members]] {:jsd jsd :cnt (count members)})))))



(def center-jsd-dist
  (future (->> ref77+PG30-SP-clusters
               #_(filter (fn[[k m]] (-> m :members count (> 1))))
               vals
               compare-clusters)))
(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg30-center-distance-data.clj"
  (prn @center-jsd-dist))

(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg30-center-distance-dist.clj"
  (prn (->> @center-jsd-dist (sort-by first)
            (mapv (fn[[jsd members]] {:jsd jsd :cnt (count members)})))))



(def center-jsd-dist
  (future (->> @ref77pg350-fut
               #_(filter (fn[[k m]] (-> m :members count (> 1))))
               vals
               compare-clusters)))
(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg350-center-distance-data.clj"
  (prn @center-jsd-dist))

(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg350-center-distance-dist.clj"
  (prn (->> @center-jsd-dist (sort-by first)
            (mapv (fn[[jsd members]] {:jsd jsd :cnt (count members)})))))

(def center-jsd-dist
  (future
    (-> "/store/data/PanClus/Stats/ref77-pg350-center-distance-data.clj"
        slurp read-string)))




(def center-merged-jsd-dist
  (future (->> ref77-pg350-with-merged-clusters
               #_(filter (fn[[k m]] (-> m :members count (> 1))))
               vals
               compare-clusters)))

(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg350-center-merged-31-distance-data.clj"
  (prn @center-merged-jsd-dist))

(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-pg350-center-merged-31-distance-dist.clj"
  (prn (->> @center-merged-jsd-dist (sort-by first)
            (mapv (fn[[jsd members]] {:jsd jsd :cnt (count members)})))))

(def center-merged-jsd-dist
  (future
    (-> "/store/data/PanClus/Stats/ref77-pg350-center-merged-distance-data.clj"
        slurp read-string)))





;;; Merged cleaned center distances
(def center-merged-cleaned-jsd-dist
  (future (->> ref77-pg350-with-merged-clusters-cleaned
               #_(filter (fn[[k m]] (-> m :members count (> 1))))
               vals
               compare-clusters)))

(io/with-out-writer
  (-> (pams/get-params :panclus-base)
      (fs/join "Stats/ref77-pg350-center-merged-distance-data.clj"))
  (prn #_@center-merged-cleaned-jsd-dist
       last-merged-jsd-dist))

(io/with-out-writer
  (-> (pams/get-params :panclus-base)
      (fs/join "Stats/ref77-pg350-center-merged-distance-dist.clj"))
  (prn (->> #_@center-merged-cleaned-jsd-dist
            last-merged-jsd-dist
            (sort-by first)
            (mapv (fn[[jsd members]] {:jsd jsd :cnt (count members)})))))




;;; Merged cleaned and leftover recluster center distances
(def center-merged-cleaned-reclu-jsd-dist
  (future (->> ref77-pg350-with-merged-clusters-cleaned-reclu
               #_(filter (fn[[k m]] (-> m :members count (> 1))))
               vals
               compare-clusters)))

(io/with-out-writer
  (-> (pams/get-params :panclus-base)
    (fs/join "Stats/ref77-pg350-center-merged-cleaned-reclu-distance-data.clj"))
  (prn @center-merged-cleaned-reclu-jsd-dist))

(io/with-out-writer
  (-> (pams/get-params :panclus-base)
    (fs/join "Stats/ref77-pg350-center-merged-cleaned-reclu-distance-dist.clj"))
  (prn (->> @center-merged-cleaned-reclu-jsd-dist (sort-by first)
            (mapv (fn[[jsd members]] {:jsd jsd :cnt (count members)})))))

(def center-merged-cleaned-reclu-jsd-dist
  (future
    (-> "/store/data/PanClus/Stats/ref77-pg350-center-merged-cleaned-reclu-distance-data.clj"
        slurp read-string)))





(->> all-singlet-center-dist (sort-by first)
     (coll/take-until #(-> % first (> 0.2)))
     count)
