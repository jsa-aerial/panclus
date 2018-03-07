








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
