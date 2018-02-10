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




"NC_003028/1442913-1443110/-1,SP_RS07540"
"NC_003028/1720276-1720479/1,SP_RS08960"

(defn get-dist-cnt [x & {:keys [AA] :or {AA true}}]
  (if (map? x)
    x
    (let [xform (if AA ntseq2aaseq identity)
          x (member->entry x)
          genomebase (pams/get-params [:biodb-info :genomes :base])
          basedir (if (= (str/substring x 0 3) "TVO")
                    (fs/join genomebase (pams/get-params
                                         [:biodb-info :genomes :tvoseq02]))
                    (fs/join genomebase (pams/get-params
                                         [:biodb-info :genomes :refseq77])))
          xsq (->> (bufiles/gen-name-seq x :basedir basedir)
                   second xform)
          xcnt (count xsq)
          xdist (p/probs 6 xsq)]
      [xdist xcnt])))

(defn get-dist [x & {:keys [AA] :or {AA true}}]
  (let [[dist cnt] (get-dist-cnt x :AA AA)]
    dist))

(defn pair-jsd [x y & {:keys [AA] :or {AA true}}]
  (let [xdist (get-dist x :AA AA)
        ydist (get-dist y :AA AA)]
    (it/jensen-shannon xdist ydist)))
#_[(it/jensen-shannon xdist ydist) (aln/align xsq ysq :match 1 :mmatch -1)]

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








;;; Transitive reduce
(def pair-map
  (->> @center-jsd-dist (sort-by first)
       (coll/take-until #(-> % first (> 0.4))) #_(take 100)
       (reduce (fn[M [scr pairs]]
                 (loop [M M
                        pairs pairs]
                   (if (empty? pairs)
                     M
                     (let [[pg1 pg2] (first pairs)]
                       (recur (assoc M
                                     pg1 (conj (M pg1 []) pg2)
                                     pg2 (conj (M pg2 []) pg1))
                              (rest pairs))))))
               {})
       #_(filter (fn[[k v]] (> (count v) 1)))
       (into {})))
(def testmap pair-map) [(count testmap) (count pair-map)]


(def testmap {"PGSP_09338" ["PGSP_02790" "PGSP_10573"], "PGSP_02790" ["PGSP_07438" "PGSP_09338"], "PGSP_05969" ["PGSP_10078" "PGSP_10253"], "PGSP_11067" ["PGSP_02794" "PGSP_08705"], "PGSP_01131" ["PGSP_06498" "PGSP_04354"], "PGSP_07228" ["PGSP_09953" "PGSP_11523"], "PGSP_05005" ["PGSP_03035" "PGSP_00060"], "PGSP_00060" ["PGSP_02303" "PGSP_05005"], "PGSP_07650" ["PGSP_00344" "PGSP_02907"], "PGSP_07193" ["PGSP_00278" "PGSP_09044"], "PGSP_07438" ["PGSP_02790" "PGSP_09237" "PGSP_10146" "PGSP_06341" "PGSP_00891"]})

(->> ["PGSP_02790"]
     (next-set #{} testmap) first
     (mapv #(->> (@ref77pg350-fut %) :members))
     (map count) m/sum)


(defn next-set
  [curset pair-map mems]
  (loop [mems mems
         curset curset
         curmap pair-map]
    (if (empty? mems)
      [curset curmap]
      (let [pgk (first mems)
            newmems (curmap pgk)
            nxtmap (dissoc curmap pgk)
            [newset newmap] (next-set (conj curset pgk) nxtmap newmems)]
        (recur (rest mems) newset newmap)))))

(defn next-id []
  (format "PGSP_%05d" (swap! PGSP-index inc)))

(defn tran-reduce
  [pair-map pgsp-start]
  (swap! PGSP-index (constantly pgsp-start))
  (reduce (fn[[M newmap] pgk]
            (if (not (newmap pgk))
              [M newmap]
              (let [[cluset newmap] (next-set #{} newmap [pgk])
                    clu-id (next-id)]
                [(assoc M clu-id cluset) newmap])))
          [{} pair-map] (keys pair-map)))


(def merged-sets
  (->> 13805 (tran-reduce {} testmap) first))



;;; Remove all pair-map (testmap) elements - result will be merged
;;; with the new clusters from tran-reduce
(def ref77pg350-minus
  (->> testmap keys (apply dissoc @ref77pg350-fut)))



(->> merged-sets vals
     (mapv (fn[s] (->> s (mapv #(->> (@ref77pg350-fut %) :members))
                      (map count) m/sum)))
     sort (coll/take-until #(> % 371)) count)

;;; count of new core in merged
(->> merged-sets vals
     (mapv (fn[s] (->> s (mapv #(->> (@ref77pg350-fut %) :members))
                      (map count) m/sum)))
     sort (filter #(= % 371)) count)

;;; Singletons count that are removed by merging
(->> pair-map keys (mapv #(->> (@ref77pg350-fut %) :members count))
     (filter #(= % 1)) count)






(defn merged-set->cluster
  [merged-set clustering & {:keys [ows] :or {ows 6}}]
  (let [[pgnm pgids] merged-set
        mems (->> pgids
                  (mapv #(->> (clustering %) :members))
                  (apply set/union))
        mdc-trips (mapv #(let [ent (member->entry %)
                               [dist cnt] (get-dist-cnt %)]
                           [% dist cnt])
                        mems)
        cnt (p/mean (mapv last mdc-trips))
        Q (it/hybrid-dictionary ows (mapv second mdc-trips))
        jsds (vfold (fn[[mem P _]] [mem (roundit (it/jensen-shannon P Q))])
                   mdc-trips)]

    {:center [pgnm cnt Q]
     :dists [Q]
     :jsds jsds
     :cnt {:sm cnt, :d 1, :m cnt}
     :members mems
     :name pgnm}))


(defn merged-sets->clusters
  [merged-sets clustering]
  (reduce (fn[C ms]
            (let [clu (merged-set->cluster ms clustering)]
              (assoc C (clu :name) clu)))
          {} merged-sets))


(def merged-clusters
  (merged-sets->clusters merged-sets))

(io/with-out-writer
  "/store/data/PanClus/Entropy/full-r77pg350-merged-clusters.clj"
  (prn merged-clusters))


;;; ref77pg350 minus all input clusters to merge plus new merged clusters
(def ref77pg350-minus
  (->> testmap keys (apply dissoc @ref77pg350-fut)))

(def ref77-pg350-with-merged-clusters
  (merge ref77pg350-minus merged-clusters))

(io/with-out-writer
  "/store/data/PanClus/Entropy/ref77-pg350-with-merged-clusters.clj"
  (prn ref77-pg350-with-merged-clusters))












;;; --- Multiple merges loses too much with transitive closure issue ---------
;;; --- So, this stuff is obsolete

(def merged-pair-map
  (->> @center-merged-jsd-dist (sort-by first)
       (coll/take-until #(-> % first (> 0.4))) #_(take 100)
       (reduce (fn[M [scr pairs]]
                 (loop [M M
                        pairs pairs]
                   (if (empty? pairs)
                     M
                     (let [[pg1 pg2] (first pairs)]
                       (recur (assoc M
                                     pg1 (conj (M pg1 []) pg2)
                                     pg2 (conj (M pg2 []) pg1))
                              (rest pairs))))))
               {})
       #_(filter (fn[[k v]] (> (count v) 1)))
       (into {})))

(def merged-sets-2
  (->> 15209 (tran-reduce merged-pair-map) first))

(->> merged-pair-map keys
     (mapv #(->> (ref77-pg350-with-merged-clusters %) :members count))
     (filter #(= % 1)) count)


(def merged-clusters-2
  (merged-sets->clusters merged-sets-2 ref77-pg350-with-merged-clusters))

(io/with-out-writer
  "/store/data/PanClus/Entropy/full-r77pg350-merged-clusters-2.clj"
  (prn merged-clusters-2))

(def ref77pg350-merged-minus
  (->> merged-pair-map keys (apply dissoc ref77-pg350-with-merged-clusters)))

(def ref77-pg350-with-merged-clusters-2
  (merge ref77pg350-merged-minus merged-clusters-2))

(io/with-out-writer
  "/store/data/PanClus/Entropy/ref77-pg350-with-merged-clusters-2.clj"
  (prn ref77-pg350-with-merged-clusters-2))
;;; --------------------------------------------------------------------------







(->> ref77-pg350-with-merged-clusters
     clus-jsd-data
     (keep (fn[m]
             (let [jsds (m :jsds)
                   >4 (->> jsds (coll/drop-until #(> % 0.35))
                           (mapv #(roundit % :places 2)))
                   mx (roundit (m :mx) :places 2)
                   mcnt (m :memcnt)
                   pgnm (m :pgnm)]
               (when (> mx 0.4)
                 {:mx mx :>4 (count >4)
                  :pgnm pgnm :memcnt mcnt}))))
     (take 2)
     #_(mapv #(vector (% :memcnt) (% :mx) (-> % :>4 count)))
     #_(mapv last)
     #_m/sum)

