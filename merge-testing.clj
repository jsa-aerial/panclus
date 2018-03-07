


;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Merging and Cleaning ;;;;;;;;;;;;;;;;;;;;;;;;;;;


(def testmap {"PGSP_09338" ["PGSP_02790" "PGSP_10573"], "PGSP_02790" ["PGSP_07438" "PGSP_09338"], "PGSP_05969" ["PGSP_10078" "PGSP_10253"], "PGSP_11067" ["PGSP_02794" "PGSP_08705"], "PGSP_01131" ["PGSP_06498" "PGSP_04354"], "PGSP_07228" ["PGSP_09953" "PGSP_11523"], "PGSP_05005" ["PGSP_03035" "PGSP_00060"], "PGSP_00060" ["PGSP_02303" "PGSP_05005"], "PGSP_07650" ["PGSP_00344" "PGSP_02907"], "PGSP_07193" ["PGSP_00278" "PGSP_09044"], "PGSP_07438" ["PGSP_02790" "PGSP_09237" "PGSP_10146" "PGSP_06341" "PGSP_00891"]})

[0.0931 [[["PGSP_06341" "PGSP_07438"] [6 4] (207.0 1305.0) 0.16]]]
[0.0759 [[["PGSP_10146" "PGSP_07438"] [3 4] (290.0 1305.0) 0.22]]]

(->> ["PGSP_06341" "PGSP_10146" "PGSP_07438"]
     (mapv #(vector % (-> (@ref77pg350-fut %) :members)))
     (mapv (fn[[pg members]]
             (vector pg (mapv #(-> % member->sq count) members)))))

(->> ["PGSP_06341" "PGSP_10146" #_"PGSP_07438"]
     (mapv #(vector % (-> (@ref77pg350-fut %) :members)))
     (mapv (fn[[pg members]] (vector pg (mapv #(-> % member->dist) members))))
     (mapcat second) (cmb/combins 2)
     (mapv (fn[[d1 d2]] (pair-jsd d1 d2)))
     sort)

(->> ["PGSP_06341" "PGSP_10146" "PGSP_07438"]
     (mapcat #(-> (@ref77pg350-fut %) :members))
     (mapv member->dist)
     (it/hybrid-dictionary 6)
     ((fn[Q]
        (->> ["PGSP_06341" "PGSP_10146" "PGSP_07438"]
             (mapv #(vector % (-> (@ref77pg350-fut %) :members)))
             (mapv (fn[[pg members]]
                     (vector pg (mapv #(-> % member->dist) members))))
             (mapcat second)
             (mapv #(pair-jsd % Q)))))
     sort)

(->> ["PGSP_02790"]
     (next-set #{} testmap) first
     (mapv #(->> (@ref77pg350-fut %) :members))
     (map count) m/sum)





(def pair-map (center-dist->pair-map jsa 0.4))

(->> #_@center-jsd-dist
     #_@center-merged-jsd-dist
     #_@center-current-merged-jsd-dist
     jsa
     (sort-by first)
     (coll/take-until #(-> % first (> 0.4))) count)


(->> pair-map keys
     (mapv #(-> (@ref77pg350-fut %) :members count))
     (group-by identity)
     (mapv (fn[[k v]] [k (count v)]))
     (sort-by first) (mapv (fn[[mcnt cnt]] (* mcnt cnt)))
     m/sum)

(->> pair-map vals (group-by #(count %)) (mapv (fn[[k v]] [k (count v)])))


(->> pair-map (group-by (fn[[k v]] (count v)))
     ((fn[m] (m 1))) (mapv (fn[[k v]] [k (first v)]))
     (mapv (fn[[c1 c2]] [[c1 (-> (@ref77pg350-fut c1) :members first)]
                        [c2 (-> (@ref77pg350-fut c2) :members first)]]))
     (vfold (fn[[[c1 m1] [c2 m2]]] [c1 c2 (roundit (pair-jsd m1 m2))]))
     (sort-by last)
     (coll/drop-until #(-> % last (> 0.15)))
     (coll/take-until #(-> % last (> 0.35)))
     count)

(->> pair-map (group-by (fn[[k v]] (count v)))
     ((fn[m] (m 1))) (mapv (fn[[k v]] [k (first v)]))
     (into {})
     (#(apply dissoc % (vals %)))
     (mapv (fn[[c1 c2]] [[c1 (-> (@ref77pg350-fut c1) :members first)]
                        [c2 (-> (@ref77pg350-fut c2) :members first)]]))
     (vfold (fn[[[c1 m1] [c2 m2]]] [c1 c2 (roundit (pair-jsd m1 m2))]))
     (sort-by last)
     #_(coll/drop-until #(-> % last (> 0.15)))
     (coll/take-until #(-> % last (> 0.35)))
     count)

(->> pair-map (group-by (fn[[k v]] (count v)))
     ((fn[m] (m 1))) (mapv (fn[[k v]] [k (first v)]))
     (into {})
     (#(apply dissoc % (vals %)))
     (mapv (fn[[c1 c2]] [[c1 (-> (@ref77pg350-fut c1) :members first)]
                        [c2 (-> (@ref77pg350-fut c2) :members first)]]))
     (vfold (fn[[[c1 m1] [c2 m2]]] [c1 c2 (roundit (pair-jsd m1 m2))]))
     (sort-by last)
     (coll/drop-until #(-> % last (> 0.15)))
     (coll/take-until #(-> % last (> 0.35)))
     (mapv second)
     (apply dissoc pair-map)
     vals (group-by #(count %)) (mapv (fn[[k v]] [k (count v)])))









;;; 13805 ref77pg350-fut version
;;; Make sure PGSP-index is correct for this!!
(def merged-sets
  (->> @PGSP-index (tran-reduce pair-map) first))

;;; count of new core in merged
(->> merged-sets vals
     (mapv (fn[s] (->> s (mapv #(->> (@ref77pg350-fut %) :members))
                      (map count) m/sum)))
     sort (filter #(= % 371)) count)

;;; Singletons count that are removed by merging
(->> pair-map keys (mapv #(->> (@ref77pg350-fut %) :members count))
     (filter #(= % 1)) count)


;;; Perform merging
;;;
(def merged-clusters
  (merged-sets->clusters merged-sets @ref77pg350-fut))

(io/with-out-writer
  "/store/data/PanClus/Entropy/r77pg350-merged-clusters-jsd-31.clj"
  (prn merged-clusters))


;;; ref77pg350 minus all input clusters to merge (remove all
;;; pair-map elements) then add back new merged clusters
(def ref77pg350-minus
  (->> pair-map keys (apply dissoc @ref77pg350-fut)))

(def ref77-pg350-with-merged-clusters
  (merge ref77pg350-minus merged-clusters))

(io/with-out-writer
  "/store/data/PanClus/Entropy/ref77-pg350-merged-jsd-31.clj"
  (prn ref77-pg350-with-merged-clusters))



(->> #_merged-clusters
     next-merged-clusters
     vals (filter #(-> % :outjsds seq))
     (mapv #(vector (-> % :jsds count) (-> % :outjsds count)))
     (sort-by second) (mapv second) m/sum)
(->> #_merged-clusters
     next-merged-clusters
     vals (filter #(-> % :outjsds seq))
     (mapv #(vector (% :name) (->> % :jsds (sort-by second) last)))
     (filter #(-> % second second (> 0.4)))
     (filter #(merged-clusters (first %))))




;;; Iterate merging... Basically just do the above steps to previous
;;; outputs repeatedly until converge.
;;;

(def current-merged-clustering ; init step
  ref77-pg350-with-merged-clusters)

#_(swap! PGSP-index
         (constantly
          (->> #_@ref77pg350-fut
               #_ref77-pg350-with-merged-clusters
               jsa-merged-clustering
               #_current-merged-clustering
               keys sort last (str/split #"_") second Integer.)))
;;; Make sure PGSP-index is correct for this!!
(def next-merged-sets
  (->> @PGSP-index (tran-reduce pair-map) first))

(def next-merged-clusters
  (merged-sets->clusters next-merged-sets current-merged-clustering))

(io/with-out-writer
  "/store/data/PanClus/Entropy/r77pg350-next-merged-clusters-jsd-39.clj"
  (prn next-merged-clusters))


;;; current clustering minus all input clusters to merge (remove all
;;; pair-map elements) then add back new merged clusters
(def current-clustering-minus
  (->> pair-map keys (apply dissoc current-merged-clustering)))

(def current-merged-clustering
  (merge current-clustering-minus next-merged-clusters))

(io/with-out-writer
  "/store/data/PanClus/Entropy/ref77-pg350-current-merged-jsd-39.clj"
  (prn current-merged-clustering))




(let [center-jsd-dist @center-current-merged-jsd-dist
      new-keys (-> next-merged-clusters keys set)
      cur-keys (-> current-merged-clustering keys set)]
  (->> center-jsd-dist (sort-by first)
       (coll/drop-until #(-> % first (> 0.6)))
       (mapcat (fn[[scr members]] (mapcat identity members)))
       (into #{})
       (set/difference cur-keys)
       #_(set/union new-keys)
       (mapv #(vector % (current-merged-clustering %)))
       count))



(def center-current-merged-jsd-dist ; current is for 0.39
  (let [center-jsd-dist @center-current-merged-jsd-dist
        new-keys (-> next-merged-clusters keys set)
        cur-keys (-> current-merged-clustering keys set)]
    (future (->> center-jsd-dist (sort-by first)
                 (coll/drop-until #(-> % first (> 0.6)))
                 (mapcat (fn[[scr members]] (mapcat identity members)))
                 (into #{})
                 (set/difference cur-keys)
                 #_(set/union new-keys)
                 (mapv #(vector % (current-merged-clustering %)))
                 (into {})
                 vals
                 compare-clusters))))

(def jsa (concat (center-distances (@DBG :current-merged-clustering)
                                   (@DBG :next-merged-clusters)
                                   @center-current-merged-jsd-dist
                                   0.6)
                 (->> @center-current-merged-jsd-dist
                      (sort-by first)
                      (coll/drop-until #(-> % first (> 0.6))))))

(def jsa-merged-clustering
  (let [pair-map (center-dist->pair-map jsa 0.4)]
    (->> pair-map
         (cur->next-merged-clusters current-merged-clustering @PGSP-index)
         (cur->next-clustering current-merged-clustering pair-map))))

(def jsa2 (concat (center-distances (@DBG :current-merged-clustering)
                                    (@DBG :next-merged-clusters)
                                    jsa
                                    0.6)
                  (->> jsa
                       (sort-by first)
                       (coll/drop-until #(-> % first (> 0.6))))))

(def jsa2-merged-clustering
  (let [pair-map (center-dist->pair-map jsa2 0.4)]
    (->> pair-map
         (cur->next-merged-clusters jsa-merged-clustering @PGSP-index)
         (cur->next-clustering jsa-merged-clustering pair-map))))



(io/with-out-writer
  "/store/data/PanClus/Stats/current-merged-39-distance-data.clj"
  (prn @center-current-merged-jsd-dist))

(io/with-out-writer
  "/store/data/PanClus/Stats/current-merged-39-distance-dist.clj"
  (prn (->> @center-current-merged-jsd-dist (sort-by first)
            (mapv (fn[[jsd members]] {:jsd jsd :cnt (count members)})))))


(+
 (->> @center-current-merged-jsd-dist (sort-by first)
      (coll/drop-until #(-> % first (> 0.6)))
      (mapcat (fn[[scr members]] (mapcat identity members)))
      (into #{}) ; count) ; 0.9 66%, 0.6 80%, 0.5 87% 0.4 95%
      (set/intersection (-> current-merged-clustering keys set))
      count)
 (->> jsa (sort-by first)
      (mapcat (fn[[scr members]] (mapcat identity members)))
      (into #{}) count))





;;; Figure out clean up post merging. The issue is JSD not satisfying
;;; triangle inequality can produce bad membership on transitive
;;; closure cluster merging.
;;;
(->> #_@ref77pg350-fut
     #_ref77pg350-minus
     #_merged-clusters
     #_next-merged-clusters
     #_ref77-pg350-with-merged-clusters
     current-merged-clustering
     clus-jsd-data
     (keep (fn[m]
             (let [x 0.4
                   jsds (m :jsds)
                   >x (->> jsds sort (coll/drop-until #(> % x))
                           (mapv #(roundit % :places 2)))
                   mx (roundit (m :mx) :places 2)
                   mcnt (m :memcnt)
                   pgnm (m :pgnm)]
               (when (> mx x)
                 {:mx mx :>x >x
                  :pgnm pgnm :memcnt mcnt})))) ; count)
     ;;(mapv :pgnm) sort count
     (mapv #(count (% :>x)))
     m/sum)


(def cleanup-ids
  (->> ref77-pg350-with-merged-clusters
     clus-jsd-data
     (keep (fn[m]
             (let [mx (roundit (m :mx) :places 2)
                   pgnm (m :pgnm)]
               (when (> mx 0.4) pgnm))))
     sort))

(def clusters-to-cleanup
  (->> cleanup-ids
       (mapv #(ref77-pg350-with-merged-clusters %))
       (mapv #(vector (% :name) %))
       (into {})))


(->> clusters-to-cleanup vals first
     ((fn[clu]
        (let [del-jsds (->> clu :jsds (sort-by second)
                            (coll/drop-until #(-> % second (> 0.35))))
              keep-jsds (->> clu :jsds (sort-by second)
                             (coll/take-until #(-> % second (> 0.35))))]
          [(assoc clu :jsds keep-jsds) del-jsds]))))



(defn new-clu-del-jsds
  [clu jsdct]
  (let [del-mems (->> clu :jsds (sort-by second)
                      (coll/drop-until #(-> % second (> jsdct)))
                      (mapv first))
        keep-jsds (->> clu :jsds (sort-by second)
                       (coll/take-until #(-> % second (> jsdct)))
                       (mapv first) (mapv #(vector % (get-dist-cnt %))))
        newcnt (p/mean (mapv #(->> % second second) keep-jsds))
        keep-jsds (mapv #(vector (first %) (-> % second first)) keep-jsds)
        centroid (it/hybrid-dictionary 6 (mapv second keep-jsds))
        newjsds (->> keep-jsds
                     (mapv #(vector (first %) (it/jensen-shannon
                                               (second %) centroid))))
        newmembers (->> newjsds (mapv first) (into #{}))
        [pgnm cnt _] (clu :center)]
    [(assoc clu
            :center [pgnm cnt centroid]
            :jsds newjsds
            :cnt {:sm newcnt, :d 1, :m newcnt}
            :dists [centroid]
            :members newmembers)
     del-mems]))

(defn clean-clusters
  [clusters jsdct]
  (vfold #(new-clu-del-jsds % jsdct)
         clusters))




(def cleaned-clus-removed-members
  (future (clean-clusters (vals clusters-to-cleanup) 0.35)))

(def cleaned-clus-removed-members @cleaned-clus-removed-members)

(io/with-out-writer
  "/store/data/PanClus/Entropy/cleaned-clusters-removed-members.clj"
  (prn cleaned-clus-removed-members))

(->> cleaned-clus-removed-members (mapcat second)
     (partition-all 2000))




(def ref77-pg350-with-merged-clusters-cleaned
  (-> (apply dissoc ref77-pg350-with-merged-clusters cleanup-ids)
      (merge (->> cleaned-clus-removed-members (mapv first)
                  (mapv #(vector (% :name) %)) (into {})))))

(def ref77-pg350-with-merged-clusters-cleaned-reclu
  (future (run-strain-clustering
           (->> cleaned-clus-removed-members (mapcat second)
                (partition-all 2000))
           ref77-pg350-with-merged-clusters-cleaned
           :pgsp-start @PGSP-index :jsdctpt 0.4 :%-max 0.60
           :chunk-size 2 :diffcut 0)))
(def ref77-pg350-with-merged-clusters-cleaned-reclu
  @ref77-pg350-with-merged-clusters-cleaned-reclu)



(swap! PGSP-index
       (fn[& args]
         (->> ref77-pg350-with-merged-clusters-cleaned
              keys sort last (str/split #"_") second Integer.)))






;;; Try clustering singlets + removed members
(->> ref77-pg350-with-merged-clusters-cleaned
     vals (filter #(-> % :members count (= 1)))
     (mapv #(-> % :members first)) (take 3))

(->> ref77-pg350-with-merged-clusters-cleaned
     vals (filter #(-> % :members count (= 1)))
     (mapv #(-> % :members first))
     (concat (mapcat second cleaned-clus-removed-members))
     panclus.clustering/gen-strain-dists)

(def jsafut nil)
  (future
    (panclus.clustering/cluster-leftovers
     (->> ref77-pg350-with-merged-clusters-cleaned
          vals (filter #(-> % :members count (= 1)))
          (mapv #(-> % :members first))
          (concat (mapcat second cleaned-clus-removed-members))
          panclus.clustering/gen-strain-dists)
     :%-max 0.7)))


(def jsafut
  (future
    (let [dists (->> ref77-pg350-with-merged-clusters-cleaned
                     vals (filter #(-> % :members count (= 1)))
                     (mapv #(-> % :members first))
                     (concat (mapcat second cleaned-clus-removed-members))
                     (partition-all 2000))
          centers (->> ref77-pg350-with-merged-clusters-cleaned
                       vals (filter #(-> % :members count (= 1)))
                       (mapv #(% :name))
                       (apply dissoc ref77-pg350-with-merged-clusters-cleaned))]
      (println :dists-chunks (count dists) :centers (count centers))
      (panclus.clustering/run-strain-clustering
       dists
       centers
       :pgsp-start @PGSP-index :jsdctpt 0.4 :%-max 0.60
       :chunk-size 2 :diffcut 0))))





