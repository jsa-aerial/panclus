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










(->> [0.0 0.0 0.0 0.08275862068965517 0.041379310344827586 0.061867002004219476 0.023694766090896352 0.042219863173742475 0.04877188529218581 0.02108631601148583 0.02108631601148583 0.00629546440276047 0.00629546440276047 0.0020328943285097666 0.043413978422945886 0.008568853176403277 0.0474951282943184 0.010548057934976188 0.010548057934976188 0.04462949697842543 0.003250120937211874]
     (mapv #(roundit % :places 2))
     (group-by identity)
     (mapv (fn[[k v]] {:x k :y (count v)})))


(->> (clus-jsd-data @ref77pg350-fut :cnt 370)
     (filter #(< (-> % :jsds count) 372))
     (random-sample 0.05)
     (map :jsds)
     (mapv (fn[x] (mapv #(roundit % :places 2) x)))
     (mapv (fn[x] (group-by identity x)))
     (mapv (fn[x] (mapv (fn[[k v]] {:x k :y (count v)}) x))))




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









(->> @center-jsd-dist first second
     (mapv (fn[pair] (mapv #(-> % ref77-SP-clusters :members first) pair)))
     (mapv (fn[[m1 m2]]
             (let [x (member->entry m1)
                   y (member->entry m2)
                   xsq (->> x bufiles/gen-name-seq second ntseq2aaseq)
                   ysq (->> y bufiles/gen-name-seq second ntseq2aaseq)
                   xdist (p/probs 6 xsq)
                   ydist (p/probs 6 ysq)]
               [m1 m2 (count xsq) (count ysq)
                (it/jensen-shannon xdist ydist)]))))


(->> center-jsd-dist (sort-by first)
     (coll/take-until #(-> % first (> 0.3)))
     (mapv (fn[[scr pairs]]
             [scr (mapv
                   (fn[pair]
                     (let [mems (mapv #(-> % ref77-SP-clusters :members)
                                      pair)
                           memcnt (mapv count mems)
                           lens (mapv #(-> % ref77-SP-clusters :center second)
                                      pair)]
                       [memcnt (m/sum memcnt)
                        lens (->> lens (apply /) double)
                        mems]))
                   pairs)]))
     clojure.pprint/pprint)

(->> center-jsd-dist (sort-by first)
     (coll/take-until #(-> % first (> 0.3)))
     (mapv (fn[[scr pairs]]
             [scr (mapv
                   (fn[pair]
                     (let [mems (mapv #(-> % ref77-SP-clusters :members)
                                      pair)
                           memcnt (mapv count mems)
                           lens (->> pair
                                     (mapv #(-> % ref77-SP-clusters
                                                :center second))
                                     sort)]
                       [memcnt (m/sum memcnt)
                        lens (->> lens (apply /) double)
                        mems]))
                   pairs)]))
     (mapv (fn[[scr pairs]] [scr (mapv (fn[[_ _ _ x]] (roundit x)) pairs)]))
     (mapv #(-> % second count))
     m/sum)


;;; All clusstering
(->> @center-jsd-dist (sort-by first)
     (coll/take-until #(-> % first (> 0.3)))
     (mapv (fn[[scr pairs]]
             [scr (mapv
                   (fn[pair]
                     (let [mems (mapv #(-> (@ref77pg350-fut %) :members)
                                      pair)
                           memcnt (mapv count mems)
                           lens (->> pair
                                     (mapv #(-> (@ref77pg350-fut %)
                                                :center second))
                                     sort)]
                       [memcnt (m/sum memcnt)
                        lens (->> lens (apply /) double)
                        mems]))
                   pairs)]))
     (mapv (fn[[scr pairs]] [scr (mapv (fn[[_ _ _ x]] (roundit x)) pairs)]))
     clojure.pprint/pprint)



(->> tigrcenter-jsd-dist (sort-by first)
     (coll/take-until #(-> % first (> 0.3)))
     (mapv (fn[[scr pairs]]
             [scr (mapv
                   (fn[pair]
                     (let [mems (mapv #(-> % tigr4-centers :members)
                                      pair)
                           memcnt (mapv count mems)
                           lens (mapv #(-> % tigr4-centers :center second)
                                      pair)]
                       [memcnt (m/sum memcnt)
                        lens (->> lens (apply /) double)
                        mems]))
                   pairs)]))
     clojure.pprint/pprint)






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












(->> ref77-SP-clusters vals
     (filter #(->> % :members (some (fn[m] (= (member->locus m) "SP_RS08960")))))
     first :members)

(->> ref77-SP-clusters vals
     (filter #(->> % :members (some (fn[m] (= (member->locus m) "SP_RS07540")))))
     first :jsds)






(->> #{#_"NC_003028/319495-319842/-1,SP_RS01685"
       #_"NC_003028/275290-275637/1,SP_RS01460"
       #_"NC_003028/938452-938799/1,SP_RS04935"
       #_"NC_003028/769180-769530/1,SP_RS04005"
       #_(->> "PGSP_11929" tigr4-centers :center last)
       ;;
       "NC_003028/1074433-1074849/-1,SP_RS05675"
       "NC_003028/1921596-1922012/-1,SP_RS10200"
       "NC_003028/20404-20820/1,SP_RS00105"
       (->> "PGSP_11343" tigr4-centers :center last)}
     (aerial.utils.math.combinatorics/combins 2)
     (mapv (fn[[x y]] (let [s1 (if (string? x) x :dist)
                           s2 (if (string? y) y :dist)]
                       [s1 s2 (pair-jsd x y)])))
     clojure.pprint/pprint)





(let [xdist (->> ref77-SP-clusters vals
                 (filter #(->> % :members
                               (some (fn[m] (= (member->locus m)
                                              "SP_RS03400")))))
                 first :center last)
      ydist (->> ref77-SP-clusters vals
                 (filter #(->> % :members
                               (some (fn[m] (= (member->locus m)
                                              "SP_RS03405")))))
                 first :center last)]
  (it/jensen-shannon xdist ydist))


(let [x (member->entry "NC_003028/769180-769530/1,SP_RS04005")
      y (member->entry "NC_003028/938452-938799/1,SP_RS04935")
      xsq (->> x bufiles/gen-name-seq second ntseq2aaseq)
      ysq (->> y bufiles/gen-name-seq second ntseq2aaseq)
      xdist (p/probs 6 xsq)
      ydist (->> test-centers make-hybrids vals
                 (filter #(->> % :members
                               (some (fn[m]
                                       (= (member->locus m) "SP_RS05675")))))
                 first :dists first)]
  (it/jensen-shannon xdist ydist))










(def tigrcenter-jsd-dist
  (->> tigr4-centers
       #_(filter (fn[[k m]] (-> m :members count (= 1))))
       vals
       panclus.stats/compare-clusters))

(->> tigr4-centers vals
     (filter #(->> % :members
                   (some (fn[m]
                           (= (member->locus m) "SP_RS05675")))))
     first :members)
#{"NC_003028/1074433-1074849/-1,SP_RS05675" "NC_003028/1921596-1922012/-1,SP_RS10200" "NC_003028/20404-20820/1,SP_RS00105"}


(->> (finish-centers test-centers 0.8 1.2) vals
     (filter #(->> % :members
                   (some (fn[m]
                           (= (member->locus m) "SP_RS05675")))))
     first :members)
#{"NC_003028/1074433-1074849/-1,SP_RS05675" "NC_003028/319495-319842/-1,SP_RS01685" "NC_003028/275290-275637/1,SP_RS01460" "NC_003028/938452-938799/1,SP_RS04935" "NC_003028/769180-769530/1,SP_RS04005"}


(->> test-centers vals
     (filter #(->> % :members
                   (some (fn[m]
                           (= (member->locus m) "SP_RS04005")))))
     first :members)
#{"NC_003028/769180-769530/1,SP_RS04005"}


(->> tigr4-centers vals
     (filter #(->> % :members
                   (some (fn[m]
                           (= (member->locus m) "SP_RS04005")))))
     first :members)
#{"NC_003028/319495-319842/-1,SP_RS01685" "NC_003028/275290-275637/1,SP_RS01460" "NC_003028/938452-938799/1,SP_RS04935" "NC_003028/769180-769530/1,SP_RS04005"}








(let [x (time (init-centers (-> reffnas first (gen-strain-dists :ows (ows :aa)))))
      cnt (->> x vals
               (filter #(-> % :members count (> 1))) count)
      jsds (->> x vals
                (filter #(-> % :members count (> 1)))
                (mapv #(-> % :jsds))
                (mapv (fn[pairs] (mapv second pairs))) flatten sort)]
  [(count x) cnt jsds])

(->> test-centers vals
     (filter #(-> % :members count (> 1)))
     (mapv #(-> % :jsds))
     (mapv (fn[pairs] (mapv second pairs))) flatten sort)
(->> test-centers vals
     (filter #(-> % :members count (> 1)))
     (mapv #(-> % :jsds))
     (mapv (fn[pairs]
             (mapv #(->> % first member->entry
                         bufiles/entry-parts second (apply -)
                         Math/abs)
                   pairs))))








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
               current-merged-clustering
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


(let [xkeys (->> (center-dist->pair-map
                  @center-current-merged-jsd-dist 0.6 :dir :above)
                 (set/intersection (-> pair-map keys set))
                 (mapcat #(conj (pair-map %) %))
                 (into #{}))
      m (center-dist->pair-map @center-current-merged-jsd-dist 0.6 :dir :above)]
  (->> xkeys (mapcat #(m %)) (into #{})
       (set/intersection (-> current-merged-clustering keys set))
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
(def jsa
  (time (let [xkeys (->> (center-dist->pair-map
                          @center-current-merged-jsd-dist 0.6 :dir :above)
                         (set/intersection (-> pair-map keys set))
                         (mapcat #(conj (pair-map %) %))
                         (into #{}))
              m (center-dist->pair-map
                 @center-current-merged-jsd-dist 0.6 :dir :above)]
          (->> xkeys (mapcat #(m %)) (into #{})
               (set/intersection (-> current-merged-clustering keys set))
               (set/union
                (let [center-jsd-dist @center-current-merged-jsd-dist
                      new-keys (-> next-merged-clusters keys set)
                      cur-keys (-> current-merged-clustering keys set)]
                  (->> center-jsd-dist (sort-by first)
                       (coll/drop-until #(-> % first (> 0.6)))
                       (mapcat (fn[[scr members]] (mapcat identity members)))
                       (into #{})
                       (set/difference cur-keys)
                       )))
               (mapv #(vector % (current-merged-clustering %)))
               (into {})
               vals
               compare-clusters))))

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
 (->> jsa2 (sort-by first)
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








;;; Find old locus (that no longer is an annotation) by coordinate search
(->> (bufiles/read-seqs
      (fs/join (pams/get-params :panclus-base) "RefSeq77/NC_003028.fna")
      :info :both)
     (map #(-> % first id-line->member))
     (map #(vector (member->locus %)
                   (->> % member->entry bufiles/entry-parts second)))
     (filter (fn[[lt [s e]]]
               (and (< (- 1678227 100) s (+ 1692710 100))
                    (< (- 1678227 100) s (+ 1692710 100))))))
=>
[["SP_RS11500" [1678222 1682844]] ["SP_RS11505" [1685971 1690224]]]

(- 1692710 1678227)
=> 14483

(map (fn[[nm [s e]]] [nm (- e s)])
     [["SP_RS11500" [1678222 1682844]] ["SP_RS11505" [1685971 1690224]]])
=>
[["SP_RS11500" 4622] ["SP_RS11505" 4253]]

(- 1685971 1682844)
=> 3127


(first (next-set #{} pair-map ["PGSP_12273"]))
(->> (next-set #{} pair-map ["PGSP_00891"])
     first (mapv #(vector % (->> (@ref77pg350-fut %) :members)))
     clojure.pprint/pprint)


(->> [["PGSP_10573" 457] ["PGSP_09338" 1019] ["PGSP_02790" 2180]
      ["PGSP_08950" 375] ["PGSP_12273" 112] ["PGSP_10146" 290]
      ["PGSP_06341" 207] ["PGSP_05682" 701] ["PGSP_09237" 593]
      ["PGSP_00891" 1518] ["PGSP_09406" 154] ["PGSP_07438" 1305]]
     (mapv (fn[[nm cnt]]
             [nm cnt (roundit
                      (pair-jsd
                       (->> (@ref77pg350-fut nm) :center last)
                       (->> (@ref77pg350-fut "PGSP_02790") :center last)))]))
     clojure.pprint/pprint)




;;; 9406 with 2790
(let [sm-9406 "TVO_1902148/2093373-2093873/1,SPS229_28830"
      pgsp-2790 ["NC_011900/1726475-1733566/1,SPN23F_RS08990"
                 "NC_012466/1593278-1599985/1,SPJ_RS08290"
                 "NC_010380/1751327-1758139/1,SPH_RS12000"
                 "TVO_1902219/1951016-1957267/1,SPS300_24225"
                 "TVO_1902008/1703275-1709589/1,SPS089_12675"
                 "TVO_1901938/1672010-1678318/1,SPS019_11040"
                 "TVO_1902191/1813537-1819851/1,SPS272_18080"]]
  (mapv #(let [[smdist smcnt] (get-dist-cnt sm-9406)
               [dist cnt] (get-dist-cnt %)
               ratio (-> cnt (/ smcnt) double)]
           [ratio (roundit (pair-jsd smdist dist))])
        pgsp-2790))


(io/with-out-writer
  (fs/join (pams/get-params :panclus-base)
           "Entropy/TVO_1902148-SPS229_28830-word-dist.clj")
  (prn (->> (member->dist "TVO_1902148/2093373-2093873/1,SPS229_28830")
            (sort-by first) vec)))

(io/with-out-writer
  (fs/join (pams/get-params :panclus-base)
           "Entropy/NC_012466-SPJ_RS08290-word-dist.clj")
  (prn (->> (member->dist "NC_012466/1593278-1599985/1,SPJ_RS08290")
            (filter #(-> % second (> 0.00099))) (sort-by first) vec)))
