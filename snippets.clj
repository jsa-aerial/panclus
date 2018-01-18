(->> "/home/jsa/Bio/ribosomal-proteins-all.csv" slurp
     csv/read-csv rest
     (mapv #(let [name (first %)
                  [s e st] (->> % (drop 2) (take 3))
                  len (Math/abs (- (Integer. e) (Integer. s)))
                  rp (last %)]
              [name len s e st rp]))
     (sort-by second >)
     (coll/takev 21))


(->> "/home/jsa/Bio/S1-proteins.csv" slurp
     csv/read-csv rest
     (mapv #(let [name (first %)
                  [s e st] (->> % (drop 2) (take 3))
                  len (Math/abs (- (Integer. e) (Integer. s)))
                  rp (last %)]
              [name s e st]))
     (sort-by first)
     (mapv #(apply bufiles/make-entry %)))

(->> "/home/jsa/Bio/S1-proteins.csv" slurp
     csv/read-csv rest
     (mapv #(let [name (first %)
                  [s e st] (->> % (drop 2) (take 3))
                  len (Math/abs (- (Integer. e) (Integer. s)))
                  rp (last %)]
              [name s e st]))
     (sort-by first)
     (mapv #(apply bufiles/make-entry %))
     bufiles/gen-name-seq-pairs
     (#(opt-wz % :alpha "AUGC"))
     #_(mapv (fn[[n s]] [n (aerial.bio.utils.seqxform/ntseq->aaseq s)]))
     #_(#(opt-wz % :alpha "ACDEFGHIKLMNPQRSSETVWY"))
     first)
;; => 10
(->> "/home/jsa/Bio/S1-proteins.csv" slurp
     csv/read-csv rest
     (mapv #(let [name (first %)
                  [s e st] (->> % (drop 2) (take 3))
                  len (Math/abs (- (Integer. e) (Integer. s)))
                  rp (last %)]
              [name s e st]))
     (sort-by first)
     (mapv #(apply bufiles/make-entry %))
     bufiles/gen-name-seq-pairs
     #_(mapv (fn[[n s]] [n (aerial.bio.utils.seqxform/ntseq->aaseq s)]))
     (aerial.utils.math.combinatorics/combins 2)
     #_shuffle
     #_(take 10)
     (mapv (fn[[nsq1 nsq2]]
             (let [[n1 sq1] nsq1
                   [n2 sq2] nsq2
                   k 10
                   d1 (p/probs k sq1)
                   d2 (p/probs k sq2)]
               [n1 n2 (it/jensen-shannon d1 d2)])))
     (mapv last)
     (apply max))
;; => 0.035


(->> jsds second (sort-by second)
     (coll/dropv-until #(-> % second (> 0.4)))
     (coll/takev-until #(-> % second (> 0.5)))
     (mapv (fn[[nm & _]] [nm (mmc2Hmap nm)])))

(mapv (fn[x] (->> x (sort-by second)
                 #_(coll/dropv-until #(-> % second (> 0.4)))
                 (coll/takev-until #(-> % second (> 0.3)))
                 (mapv (fn[[nm & _]] [nm (mmc2Hmap nm)])) count))
      jsds)

#_(->> :????
     (filter (fn[[k v]] (not= 3 (count v))))
     (filter (fn[[nm recs]] (not (some #(->> % second (= "NC_003028")) recs))))
     (mapv first))

;;; "NC_018630"
(def Strains ["NC_003028" "NC_012469" "NC_008533"])
(def mmc2 (->> "/home/jsa/Bio/mmc2.csv" slurp csv/read-csv rest))
(def mmc2H (->> mmc2 rest
                (filter (fn[[t4 t19f d39 & _]]
                          (reduce (fn[all x] (and all (not= x "n.p.")))
                                  true [t4 t19f d39])))))
(def mmc2Hmap (reduce (fn[M [t4 t19 d39]]
                        (assoc M t4 [t19 d39]
                               t19 [t4 d39]
                               d39 [t4 t19]))
                      {} mmc2H))
(def mmc2T4Hmap (reduce (fn[M [t4 t19 d39]]
                        (assoc M t4 t4
                               t19 t4
                               d39 t4))
                      {} mmc2H))

(def homologs
  (->> Strains ;; only works for t4, 19f, d39 - only xref table
       (mapcat (fn[nc]
                 (let [base "/store/data/SPStrains/RefSeq77/CSV/"]
                   (->> (str nc ".csv") (fs/join base)  slurp
                        csv/read-csv rest (mapv #(coll/takev 6 %))
                        (filter (fn[[nm s e st g ol]]
                                  (or (not= g "NA") ; gene or known xref locus
                                      (mmc2T4Hmap ol))))))))
       (map (fn[[nm s e st g ol]]
              (let [len (Math/abs (- (Integer. e) (Integer. s)))]
                [len nm s e st g ol])))
       (group-by (fn[[len nm s e st g ol]]
                   (if (not= g "NA") g (mmc2T4Hmap ol))))
       (filter (fn[[k v]] (= 3 (count v)))) #_(take 10)))

(def OWS
  [(->> homologs (sort-by (fn[[k recs]] (ffirst recs))) (random-sample 0.05)
        (mapv (fn[[k recs]]
                [k (ffirst recs)
                 (->> recs (map #(->> % rest (take 4)
                                      (apply bufiles/make-entry)))
                      bufiles/gen-name-seq-pairs
                      (#(opt-wz % :alpha "AUGC"))
                      first)])))
   (->> homologs (sort-by (fn[[k recs]] (ffirst recs))) (random-sample 0.05)
        (mapv (fn[[k recs]]
                [k (ffirst recs)
                 (->> recs (map #(->> % rest (take 4)
                                      (apply bufiles/make-entry)))
                      bufiles/gen-name-seq-pairs
                      (mapv (fn[[n s]]
                              [n (sqx/ntseq->aaseq s)]))
                      (#(opt-wz % :alpha "ACDEFGHIKLMNPQRSSETVWY"))
                      first)])))])
;;=> NT: [["gyrB" 10] ["cbiO" 9]]
;;=> AA: [["gyrB" 6] ["cbiO" 5]]


(def jsds
  (let [samps (->> homologs (sort-by (fn[[k recs]] (ffirst recs)))
                   #_(random-sample 0.20))]
    [(->> samps
          (vfold (fn[[k recs]]
                   (let [x (->> recs (map #(->> % rest (take 4)
                                                (apply bufiles/make-entry)))
                                bufiles/gen-name-seq-pairs
                                (cmb/combins 2)
                                (mapv (fn[[nsq1 nsq2]]
                                        (let [[n1 sq1] nsq1
                                              [n2 sq2] nsq2
                                              k 9
                                              d1 (p/probs k sq1)
                                              d2 (p/probs k sq2)]
                                          [n1 n2 (it/jensen-shannon d1 d2)])))
                                (sort-by last))]
                     [k x])))
          (sort-by (comp last first second) >))
     (->> samps
          (vfold (fn[[k recs]]
                   (let [x (->> recs (map #(->> % rest (take 4)
                                                (apply bufiles/make-entry)))
                                bufiles/gen-name-seq-pairs
                                (mapv (fn[[n s]]
                                        [n (ntseq2aaseq s)]))
                                (cmb/combins 2)
                                (mapv (fn[[nsq1 nsq2]]
                                        (let [[n1 sq1] nsq1
                                              [n2 sq2] nsq2
                                              k 6
                                              d1 (p/probs k sq1)
                                              d2 (p/probs k sq2)]
                                          [n1 n2 (it/jensen-shannon d1 d2)])))
                                (sort-by last))]
                     [k x])))
          (sort-by (comp last first second) >))]))


;;; Just the _named_ (well known) genes
(->> jsds second (sort-by (comp last first second))
     (coll/takev-until #(-> % ((comp last first second)) (> 0.42)))
     (filter (fn[[n & _]] (not= "SP" (str/take 2 n)))) count)



(def ntjsds-m
  (->> jsds first
       (vfold (fn[[n v]]
                (let [v (sort-by last > v)
                      sq-pairs
                      (mapv (fn[x] (->> x butlast
                                       (mapv #(->> % bufiles/gen-name-seq
                                                   second))))
                            v)
                      hms-lens (mapv (fn[[s1 s2]]
                                       (let [hm (it/levenshtein s1 s2)
                                             len (max (count s1) (count s2))]
                                         [(-> hm (/ len) double) [hm len]]))
                                     sq-pairs)
                      v2 (->> v (interleave hms-lens)
                              (partition-all 2)
                              (mapv (fn[[x y]] (conj x y))))]
                  [n v2])))
       (into {})))

(def aajsds-m
  (->> jsds second
       (vfold (fn[[n v]]
                (let [v (sort-by last > v)
                      sq-pairs
                      (mapv (fn[x] (->> x butlast
                                       (mapv #(->> % bufiles/gen-name-seq
                                                   second))
                                       (mapv ntseq2aaseq)))
                            v)
                      hms-lens (mapv (fn[[s1 s2]]
                                       (let [hm (it/levenshtein s1 s2)
                                             len (max (count s1) (count s2))]
                                         [(-> hm (/ len) double) [hm len]]))
                                     sq-pairs)
                      v2 (->> v (interleave hms-lens)
                              (partition-all 2)
                              (mapv (fn[[x y]] (conj x y))))]
                  [n v2])))
       (into {})))


(aajsds-m "engB")
(ntjsds-m "recX")

(->> ["NC_012469/1519593-1522385/-1" "NC_018630/1506925-1509717/-1"]
     bufiles/gen-name-seq-pairs (map second) (map count))

(->> ["NC_003028/1475387-1475974/-1" "NC_012469/1439232-1439819/-1" ]
     bufiles/gen-name-seq-pairs (map second) (map count))

(->> jsds second (sort-by (comp first second))
     (coll/takev-until #(-> % ((comp first second)) (> 0.42)))
     #_(filter (fn[[n & _]] (not= "SP" (str/take 2 n))))
     (filter (fn[[_ [mn mx av] & _]] (> mx 0.99))) last)


(defn roundit [r & {:keys [places] :or {places 4}}]
  (let [n (Math/pow 10.0 places)]
    (-> r (* n) Math/round (/ n))))



(def gene2jsds
  (->> ntjsds-m keys
       (mapv #(vector % {:nt (ntjsds-m %) :aa (aajsds-m %)}))
       (into {})))


(def nm-xref {"NC_003028" "T4", "NC_012469" "19F", "NC_008533" "D39"})

(def vgl-maps
  (coll/concatv
   (->> gene2jsds
        (mapcat (fn[[k m]]
                  (map (fn[[hm% [hm len] [e1 e2 jsd]]]
                         (let [n1 (->> e1 (str/split #"/") first nm-xref)
                               n2 (->> e2 (str/split #"/") first nm-xref)
                               hm% (roundit hm%)
                               jsd (roundit jsd)]
                           {:name k :hm hm% :len len
                            :kind "NT" :nt jsd
                            :tt (format "[%s %s] %s: %s,%s %s"
                                        hm% jsd k n1 n2 [hm len])}))
                       (m :nt)))))
   (->> gene2jsds
        (mapcat (fn[[k m]]
                  (map (fn[[hm% [hm len] [e1 e2 jsd]]]
                         (let [n1 (->> e1 (str/split #"/") first nm-xref)
                               n2 (->> e2 (str/split #"/") first nm-xref)
                               hm% (roundit hm%)
                               jsd (roundit jsd)]
                           {:name k :hm hm% :len len
                            :kind "AA" :aa jsd
                            :tt (format "[%s %s] %s: %s,%s %s"
                                        hm% jsd k n1 n2 [hm len])}))
                       (m :aa)))))))


(io/with-out-writer "/home/jsa/Bio/panclus-vgl-maps.clj"
  (prn vgl-maps))






(->> ref-dists-all
     (concat (gen-strain-group-dists (pgfnas :pg30) :ows (ows :aa)))
     (apply concat) count)
110339


(io/with-out-writer "/store/data/PanClus/Entropy/ref77-SP-clusters.clj"
  (prn ref77-SP-clusters))

(->> ref77-SP-clusters vals
     (group-by #(count (% :members)))
     (map (fn[[k m]] [k (count m)]))
     (sort-by first)
     (mapv (fn[[sz ct]] {:sz sz :cnt ct})))
=>




(io/with-out-writer "/store/data/PanClus/Entropy/ref77+PG30-SP-clusters.clj"
  (prn ref77+PG30-SP-clusters))

(->> ref77+PG30-SP-clusters vals
     (group-by #(count (% :members)))
     (map (fn[[k m]] [k (count m)]))
     (sort-by first)
     (mapv (fn[[sz ct]] {:sz sz :cnt ct})))
=>
[{:sz 1, :cnt 2183} {:sz 2, :cnt 1109} {:sz 3, :cnt 444} {:sz 4, :cnt 360} {:sz 5, :cnt 309} {:sz 6, :cnt 185} {:sz 7, :cnt 168} {:sz 8, :cnt 145} {:sz 9, :cnt 125} {:sz 10, :cnt 93} {:sz 11, :cnt 54} {:sz 12, :cnt 66} {:sz 13, :cnt 58} {:sz 14, :cnt 56} {:sz 15, :cnt 41} {:sz 16, :cnt 46} {:sz 17, :cnt 43} {:sz 18, :cnt 34} {:sz 19, :cnt 47} {:sz 20, :cnt 50} {:sz 21, :cnt 34} {:sz 22, :cnt 24} {:sz 23, :cnt 23} {:sz 24, :cnt 31} {:sz 25, :cnt 26} {:sz 26, :cnt 24} {:sz 27, :cnt 20} {:sz 28, :cnt 34} {:sz 29, :cnt 24} {:sz 30, :cnt 23} {:sz 31, :cnt 24} {:sz 32, :cnt 23} {:sz 33, :cnt 21} {:sz 34, :cnt 17} {:sz 35, :cnt 16} {:sz 36, :cnt 20} {:sz 37, :cnt 24} {:sz 38, :cnt 17} {:sz 39, :cnt 16} {:sz 40, :cnt 15} {:sz 41, :cnt 16} {:sz 42, :cnt 21} {:sz 43, :cnt 21} {:sz 44, :cnt 20} {:sz 45, :cnt 29} {:sz 46, :cnt 42} {:sz 47, :cnt 33} {:sz 48, :cnt 46} {:sz 49, :cnt 74} {:sz 50, :cnt 263} {:sz 51, :cnt 915} {:sz 52, :cnt 1} {:sz 53, :cnt 1} {:sz 54, :cnt 3} {:sz 55, :cnt 1} {:sz 56, :cnt 1} {:sz 57, :cnt 1} {:sz 59, :cnt 1} {:sz 60, :cnt 2} {:sz 62, :cnt 1} {:sz 65, :cnt 1} {:sz 74, :cnt 1} {:sz 76, :cnt 1} {:sz 77, :cnt 2} {:sz 109, :cnt 2} {:sz 132, :cnt 1}]




(io/with-out-writer "/store/data/PanClus/Entropy/ref77+PG350-SP-clusters.clj"
  (prn @ref77pg350-fut))

(->> @ref77pg350-fut vals
     (group-by #(count (% :members)))
     (map (fn[[k m]] [k (count m)]))
     (sort-by first)
     (mapv (fn[[sz ct]] {:sz sz :cnt ct})))




(->> ref77-SP-clusters
     clus-jsd-data
     (group-by #(-> % :mean (roundit :places 2)))
     (map (fn[[k m]] [k (count m)]))
     (sort-by first)
     (mapv (fn[[jsd ct]] {:jsd jsd :cnt ct})))


(->> ref77+PG30-SP-clusters
     clus-jsd-data
     (group-by #(-> % :mean (roundit :places 2)))
     (map (fn[[k m]] [k (count m)]))
     (sort-by first)
     (mapv (fn[[jsd ct]] {:jsd jsd :cnt ct})))


(->> @ref77pg350-fut
     clus-jsd-data
     (group-by #(-> % :mean (roundit :places 2)))
     (map (fn[[k m]] [k (count m)]))
     (sort-by first)
     (mapv (fn[[jsd ct]] {:jsd jsd :cnt ct})))




(->> ref77+PG30-SP-clusters
     (#(clus-jsd-data % :cnt 1))
     (group-by #(-> % :mx (roundit :places 2)))
     (map (fn[[k m]] [k (count m)]))
     (sort-by first)
     (mapv (fn[[mx ct]] {:mx mx :cnt ct}))
     (filter #(<= (% :mx) 0.4))
     (mapv :cnt)
     m/sum)

(->> @ref77pg350-fut
     clus-jsd-data
     (group-by #(-> % :mx (roundit :places 2)))
     (map (fn[[k m]] [k (count m)]))
     (sort-by first)
     (mapv (fn[[mx ct]] {:mx mx :cnt ct}))
     (filter #(<= (% :mx) 0.4))
     (mapv :cnt)
     m/sum)





(clus-counts ref77-SP-clusters :sdev-grp :less+1sdev :cut% 0.8)
(clus-counts ref77-SP-clusters :sdev-grp :within-1sdev :cut% 0.8)

(clus-counts ref77+PG30-SP-clusters :sdev-grp :less+1sdev :cut% 0.8)
(clus-counts ref77+PG30-SP-clusters :sdev-grp :within-1sdev :cut% 0.8)

(clus-counts @ref77pg350-fut :sdev-grp :less+1sdev :cut% 0.8)
(clus-counts @ref77pg350-fut :sdev-grp :within-1sdev :cut% 0.8)



(->> [0.0 0.0 0.0 0.08275862068965517 0.041379310344827586 0.061867002004219476 0.023694766090896352 0.042219863173742475 0.04877188529218581 0.02108631601148583 0.02108631601148583 0.00629546440276047 0.00629546440276047 0.0020328943285097666 0.043413978422945886 0.008568853176403277 0.0474951282943184 0.010548057934976188 0.010548057934976188 0.04462949697842543 0.003250120937211874]
     (mapv #(roundit % :places 2))
     (group-by identity)
     (mapv (fn[[k v]] {:x k :y (count v)})))



(->> (clus-jsd-data @ref77pg350-fut :cnt 370)
     (filter #(< (-> % :jsds count) 372))
     (mapv #(vector (% :mean) (% :sdev)))
     (sort-by second))

(->> (clus-jsd-data @ref77pg350-fut :cnt 10)
     (filter #(< (-> % :jsds count) 372))
     (mapv #(vector (% :mean) (% :sdev)))
     (sort-by second)
     (coll/drop-until #(-> % second (> 0.1)))
     count)

(io/with-out-writer
  "/store/data/PanClus/Stats/cluster-jsd-means-sdevs-371-only.clj"
  (prn (->> (clus-jsd-data @ref77pg350-fut :cnt 370)
            (filter #(< (-> % :jsds count) 372))
            (mapv #(vector (-> % :mean roundit) (-> % :sdev roundit)))
            (sort-by second)
            #_(coll/drop-until #(-> % second (> 0.1)))
            (group-by #(-> % second))
            (mapv (fn[[sdev v]] {:sdev sdev :cnt (count v)})))))

(io/with-out-writer
  "/store/data/PanClus/Stats/cluster-jsd-means-sdevs-372-plus.clj"
  (prn (->> (clus-jsd-data @ref77pg350-fut :cnt 371)
            #_(filter #(< (-> % :jsds count) 372))
            (mapv #(vector (-> % :mean roundit) (-> % :sdev roundit)))
            (sort-by second)
            #_(coll/drop-until #(-> % second (> 0.1)))
            (group-by #(-> % second))
            (mapv (fn[[sdev v]] {:sdev sdev :cnt (count v)})))))


(->> (range 0.1 0.9 0.1) (map roundit)
     (mapv #(vector % (clus-counts @ref77pg350-fut
                                   :sdev-grp :within-1sdev :cut% %)))
     (mapv (fn[[per [cnt cper]]] {:cut per :cnt cnt :per cper})))

(->> (range 0.1 0.9 0.1) (map roundit)
     (mapv #(vector % (clus-counts @ref77pg350-fut
                                   :sdev-grp :less+1sdev :cut% %)))
     (mapv (fn[[per [cnt cper]]] {:cut per :cnt cnt :per cper})))




(->> (clus-jsd-data @ref77pg350-fut :cnt 370)
     (filter #(< (-> % :jsds count) 372))
     (random-sample 0.05)
     (map :jsds)
     (mapv (fn[x] (mapv #(roundit % :places 2) x)))
     (mapv (fn[x] (group-by identity x)))
     (mapv (fn[x] (mapv (fn[[k v]] {:x k :y (count v)}) x))))




(defn singlet-clu->distpt [clu]
  (let [nm (-> clu :members first)
        [_ cnt dist] (clu :center)]
    [nm cnt dist]))

(defn next-point-dist [clusters mode]
  (if (= mode :center-center)
    (-> clusters first :center last)
    (assert :NYI "only :center-center supported at this point")))

(defn compare-clusters [clusters & {:keys [mode] :or {mode :center-center}}]
  (loop [comp-dist {}
         clusters clusters]
    (if (empty? (rest clusters))
      comp-dist
      (let [point-dist (next-point-dist clusters mode)
            dists (mapv  #(vector (% :name)(-> % :center last)) (rest clusters))
            [nm score] (->> dists
                            (vfold (fn[[nm Q]]
                                     [nm (roundit (it/jensen-shannon point-dist Q))]))
                            (sort-by second) first)]
        (println (format "Point %s : [%s %s]" ((first clusters) :name)
                         nm score))
        (recur (assoc comp-dist score (inc (get comp-dist score 0)))
               (rest clusters))))))


(def center-jsd-dist
  (->> ref77-SP-clusters
       (filter (fn[[k m]] (-> m :members count (= 1))))
       vals
       compare-clusters))
;;; WTF?!?!
Point PGSP_05505 : [PGSP_05707 0.2048]
(->> "PGSP_05505" ref77-SP-clusters ((fn[m] [(-> m :center second) (-> m :members first)])))
(->> "PGSP_05707" ref77-SP-clusters ((fn[m] [(-> m :center second) (-> m :members first)])))

Point PGSP_05384 : [PGSP_05553 0.1655]
(->> "PGSP_05384" ref77-SP-clusters ((fn[m] [(-> m :center second) (-> m :members first)])))
(->> "PGSP_05553" ref77-SP-clusters ((fn[m] [(-> m :center second) (-> m :members first)])))

Point PGSP_05744 : [PGSP_06906 0.2457]
Point PGSP_06369 : [PGSP_06286 0.2643]

(io/with-out-writer
  "/store/data/PanClus/Stats/ref77-center-distance-dist.clj"
  (prn (->> center-jsd-dist (sort-by first)
            (mapv (fn[[jsd cnt]] {:jsd jsd :cnt cnt})))))
(io/with-out-writer
  "/store/data/PanClus/Stats/ref77+PG30-center-distance-dist.clj"
  (prn (->> center-jsd-dist (sort-by first)
            (mapv (fn[[jsd cnt]] {:jsd jsd :cnt cnt})))))
(io/with-out-writer
  "/store/data/PanClus/Stats/ref77+PG350-singlet-center-distance-dist.clj"
  (prn (->> @center-jsd-dist-fut (sort-by first)
            (mapv (fn[[jsd cnt]] {:jsd jsd :cnt cnt})))))


(def center-jsd-dist-fut
  (future (->> @ref77pg350-fut
               (filter (fn[[k m]] (-> m :members count (> 1))))
               vals
               compare-clusters)))








(defn xfn [id]
  (->> id (str/split #",")
       first bufiles/gen-name-seq
       ((fn[[n s]] [id (->> s ntseq2aaseq (p/probs 6))]))))

(defn yfn [clu]
  [(first clu) (->> clu second :members (vfold xfn) (into {}))])

(->> ref77-SP-clusters (mapv yfn))



(let [cdist (->> "PGSP_00689" ref77-SP-clusters :center last)
      [id mdists] (->> ref77-SP-clusters (take 3) (mapv yfn) first)]
  [id (sort-by second (vfold (fn[[n Q]]
                               [n (it/jensen-shannon cdist Q)])
                             mdists))])
