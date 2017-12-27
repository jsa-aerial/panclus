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





(->> ref-dists (concat (gen-strain-group-dists (pgfnas :pg30) :ows (ows :aa)))
     (apply concat) count)
108123

(->> ref77+PG30-SP-clusters vals
     #_(filter #(> (count (% :members)) 1))
     (group-by #(count (% :members)))
     (map (fn[[k m]] [k (count m)]))
     (sort-by first)
     (mapv (fn[[sz ct]] {:sz sz :cnt ct})))
=>
[{:sz 2, :cnt 26} {:sz 3, :cnt 19} {:sz 4, :cnt 18} {:sz 5, :cnt 10} {:sz 6, :cnt 24} {:sz 7, :cnt 15} {:sz 8, :cnt 18} {:sz 9, :cnt 20} {:sz 10, :cnt 7} {:sz 11, :cnt 12} {:sz 12, :cnt 7} {:sz 13, :cnt 8} {:sz 14, :cnt 10} {:sz 15, :cnt 19} {:sz 16, :cnt 12} {:sz 17, :cnt 11} {:sz 18, :cnt 26} {:sz 19, :cnt 33} {:sz 20, :cnt 23} {:sz 21, :cnt 7} {:sz 22, :cnt 11} {:sz 23, :cnt 13} {:sz 24, :cnt 13} {:sz 25, :cnt 12} {:sz 26, :cnt 11} {:sz 27, :cnt 12} {:sz 28, :cnt 14} {:sz 29, :cnt 15} {:sz 30, :cnt 11} {:sz 31, :cnt 23} {:sz 32, :cnt 13} {:sz 33, :cnt 17} {:sz 34, :cnt 15} {:sz 35, :cnt 9} {:sz 36, :cnt 12} {:sz 37, :cnt 15} {:sz 38, :cnt 11} {:sz 39, :cnt 20} {:sz 40, :cnt 10} {:sz 41, :cnt 19} {:sz 42, :cnt 20} {:sz 43, :cnt 19} {:sz 44, :cnt 28} {:sz 45, :cnt 40} {:sz 46, :cnt 31} {:sz 47, :cnt 50} {:sz 48, :cnt 79} {:sz 49, :cnt 268} {:sz 50, :cnt 935} {:sz 51, :cnt 2} {:sz 52, :cnt 2} {:sz 53, :cnt 1} {:sz 55, :cnt 1} {:sz 56, :cnt 1} {:sz 57, :cnt 1} {:sz 58, :cnt 2} {:sz 59, :cnt 4} {:sz 60, :cnt 3} {:sz 61, :cnt 2} {:sz 63, :cnt 2} {:sz 67, :cnt 4} {:sz 68, :cnt 1} {:sz 71, :cnt 4} {:sz 73, :cnt 6} {:sz 75, :cnt 1} {:sz 76, :cnt 2} {:sz 77, :cnt 1} {:sz 78, :cnt 2} {:sz 83, :cnt 2} {:sz 90, :cnt 1} {:sz 91, :cnt 2} {:sz 93, :cnt 1} {:sz 94, :cnt 1} {:sz 98, :cnt 11} {:sz 106, :cnt 5} {:sz 107, :cnt 1} {:sz 117, :cnt 2} {:sz 126, :cnt 1} {:sz 127, :cnt 1} {:sz 131, :cnt 1} {:sz 133, :cnt 1} {:sz 136, :cnt 4} {:sz 159, :cnt 5} {:sz 175, :cnt 1} {:sz 181, :cnt 1} {:sz 194, :cnt 1} {:sz 202, :cnt 1} {:sz 219, :cnt 1} {:sz 223, :cnt 1} {:sz 224, :cnt 4} {:sz 225, :cnt 2} {:sz 313, :cnt 1}]


(->> jsa (coll/concatv ref77+PG30-SP-clusters)
     (group-by #(count (% :members)))
     (map (fn[[k m]] [k (count m)]))
     (sort-by first)
     (map (fn[[sz ct]] {:sz sz :cnt ct})))
=>
[{:sz 2, :cnt 1085} {:sz 3, :cnt 390} {:sz 4, :cnt 323} {:sz 5, :cnt 282} {:sz 6, :cnt 192} {:sz 7, :cnt 140} {:sz 8, :cnt 135} {:sz 9, :cnt 122} {:sz 10, :cnt 82} {:sz 11, :cnt 55} {:sz 12, :cnt 58} {:sz 13, :cnt 46} {:sz 14, :cnt 48} {:sz 15, :cnt 41} {:sz 16, :cnt 38} {:sz 17, :cnt 33} {:sz 18, :cnt 47} {:sz 19, :cnt 50} {:sz 20, :cnt 36} {:sz 21, :cnt 21} {:sz 22, :cnt 26} {:sz 23, :cnt 22} {:sz 24, :cnt 25} {:sz 25, :cnt 24} {:sz 26, :cnt 16} {:sz 27, :cnt 22} {:sz 28, :cnt 31} {:sz 29, :cnt 20} {:sz 30, :cnt 21} {:sz 31, :cnt 31} {:sz 32, :cnt 18} {:sz 33, :cnt 25} {:sz 34, :cnt 18} {:sz 35, :cnt 12} {:sz 36, :cnt 16} {:sz 37, :cnt 25} {:sz 38, :cnt 14} {:sz 39, :cnt 23} {:sz 40, :cnt 14} {:sz 41, :cnt 21} {:sz 42, :cnt 22} {:sz 43, :cnt 22} {:sz 44, :cnt 28} {:sz 45, :cnt 42} {:sz 46, :cnt 32} {:sz 47, :cnt 50} {:sz 48, :cnt 82} {:sz 49, :cnt 269} {:sz 50, :cnt 935} {:sz 51, :cnt 3} {:sz 52, :cnt 2} {:sz 53, :cnt 1} {:sz 54, :cnt 1} {:sz 55, :cnt 1} {:sz 56, :cnt 1} {:sz 57, :cnt 1} {:sz 58, :cnt 2} {:sz 59, :cnt 4} {:sz 60, :cnt 4} {:sz 61, :cnt 2} {:sz 63, :cnt 2} {:sz 67, :cnt 4} {:sz 68, :cnt 1} {:sz 71, :cnt 4} {:sz 73, :cnt 6} {:sz 75, :cnt 1} {:sz 76, :cnt 2} {:sz 77, :cnt 1} {:sz 78, :cnt 2} {:sz 83, :cnt 2} {:sz 90, :cnt 1} {:sz 91, :cnt 2} {:sz 93, :cnt 1} {:sz 94, :cnt 1} {:sz 98, :cnt 11} {:sz 106, :cnt 5} {:sz 107, :cnt 1} {:sz 117, :cnt 2} {:sz 126, :cnt 1} {:sz 127, :cnt 1} {:sz 131, :cnt 1} {:sz 133, :cnt 1} {:sz 136, :cnt 4} {:sz 159, :cnt 5} {:sz 175, :cnt 1} {:sz 181, :cnt 1} {:sz 194, :cnt 1} {:sz 202, :cnt 1} {:sz 219, :cnt 1} {:sz 223, :cnt 1} {:sz 224, :cnt 4} {:sz 225, :cnt 2} {:sz 313, :cnt 1}]



(io/with-out-writer "/store/data/PanClus/Entropy/ref77+PG30-SP-clusters.clj"
  (prn ref77+PG30-SP-clusters))

(->> ref77+PG30-SP-clusters vals
     (group-by #(count (% :members)))
     (map (fn[[k m]] [k (count m)]))
     (sort-by first)
     (mapv (fn[[sz ct]] {:sz sz :cnt ct})))
=>
[{:sz 2, :cnt 94} {:sz 3, :cnt 101} {:sz 4, :cnt 78} {:sz 5, :cnt 124} {:sz 6, :cnt 156} {:sz 7, :cnt 121} {:sz 8, :cnt 114} {:sz 9, :cnt 170} {:sz 10, :cnt 473} {:sz 11, :cnt 1492} {:sz 12, :cnt 1102} {:sz 13, :cnt 6} {:sz 14, :cnt 18} {:sz 15, :cnt 7} {:sz 16, :cnt 5} {:sz 17, :cnt 4} {:sz 18, :cnt 12} {:sz 19, :cnt 6} {:sz 20, :cnt 3} {:sz 21, :cnt 10} {:sz 22, :cnt 10} {:sz 23, :cnt 4} {:sz 24, :cnt 5} {:sz 25, :cnt 3} {:sz 27, :cnt 8} {:sz 28, :cnt 4} {:sz 29, :cnt 6} {:sz 31, :cnt 3} {:sz 32, :cnt 13} {:sz 33, :cnt 11} {:sz 34, :cnt 1} {:sz 35, :cnt 2} {:sz 36, :cnt 2} {:sz 43, :cnt 1} {:sz 44, :cnt 8} {:sz 45, :cnt 4} {:sz 46, :cnt 6} {:sz 47, :cnt 1} {:sz 48, :cnt 1} {:sz 49, :cnt 1} {:sz 50, :cnt 2} {:sz 51, :cnt 1} {:sz 52, :cnt 1} {:sz 53, :cnt 2} {:sz 54, :cnt 1} {:sz 55, :cnt 1} {:sz 56, :cnt 1} {:sz 58, :cnt 1}]

(io/with-out-writer
  "/store/data/PanClus/Entropy/ref77+PG30-SP-clusters-centers-jsdcpt4.clj"
  (prn ref77+PG30-SP-clusters))


[{:sz 0, :cnt 35} {:sz 1, :cnt 26} {:sz 2, :cnt 134} {:sz 3, :cnt 104} {:sz 4, :cnt 108} {:sz 5, :cnt 120} {:sz 6, :cnt 146} {:sz 7, :cnt 142} {:sz 8, :cnt 133} {:sz 9, :cnt 165} {:sz 10, :cnt 439} {:sz 11, :cnt 1203} {:sz 12, :cnt 64} {:sz 13, :cnt 77} {:sz 14, :cnt 65} {:sz 15, :cnt 82} {:sz 16, :cnt 73} {:sz 17, :cnt 79} {:sz 18, :cnt 109} {:sz 19, :cnt 139} {:sz 20, :cnt 313} {:sz 21, :cnt 1038} {:sz 22, :cnt 35} {:sz 23, :cnt 48} {:sz 24, :cnt 31} {:sz 25, :cnt 42} {:sz 26, :cnt 29} {:sz 27, :cnt 39} {:sz 28, :cnt 38} {:sz 29, :cnt 29} {:sz 30, :cnt 32} {:sz 31, :cnt 39} {:sz 32, :cnt 49} {:sz 33, :cnt 34} {:sz 34, :cnt 38} {:sz 35, :cnt 42} {:sz 36, :cnt 51} {:sz 37, :cnt 64} {:sz 38, :cnt 69} {:sz 39, :cnt 101} {:sz 40, :cnt 276} {:sz 41, :cnt 964} {:sz 42, :cnt 25} {:sz 43, :cnt 27} {:sz 44, :cnt 43} {:sz 45, :cnt 48} {:sz 46, :cnt 35} {:sz 47, :cnt 48} {:sz 48, :cnt 91} {:sz 49, :cnt 270} {:sz 50, :cnt 920} {:sz 51, :cnt 14} {:sz 52, :cnt 2} {:sz 53, :cnt 7} {:sz 54, :cnt 1} {:sz 55, :cnt 3} {:sz 56, :cnt 1} {:sz 57, :cnt 5} {:sz 58, :cnt 3} {:sz 59, :cnt 5} {:sz 60, :cnt 1} {:sz 61, :cnt 1} {:sz 62, :cnt 7} {:sz 63, :cnt 2} {:sz 66, :cnt 6} {:sz 67, :cnt 2} {:sz 68, :cnt 1} {:sz 69, :cnt 2} {:sz 70, :cnt 1} {:sz 71, :cnt 2} {:sz 73, :cnt 5} {:sz 74, :cnt 1} {:sz 75, :cnt 4} {:sz 76, :cnt 4} {:sz 77, :cnt 2} {:sz 79, :cnt 2} {:sz 80, :cnt 1} {:sz 81, :cnt 1} {:sz 82, :cnt 1} {:sz 84, :cnt 1} {:sz 87, :cnt 2} {:sz 88, :cnt 2} {:sz 89, :cnt 5} {:sz 92, :cnt 2} {:sz 93, :cnt 2} {:sz 94, :cnt 2} {:sz 95, :cnt 1} {:sz 96, :cnt 1} {:sz 97, :cnt 3} {:sz 98, :cnt 12} {:sz 103, :cnt 1} {:sz 108, :cnt 5} {:sz 112, :cnt 5} {:sz 115, :cnt 2} {:sz 121, :cnt 1} {:sz 122, :cnt 1} {:sz 124, :cnt 1} {:sz 127, :cnt 1} {:sz 129, :cnt 1} {:sz 132, :cnt 2} {:sz 134, :cnt 1} {:sz 138, :cnt 1} {:sz 152, :cnt 1} {:sz 154, :cnt 1} {:sz 156, :cnt 1} {:sz 157, :cnt 1} {:sz 158, :cnt 1} {:sz 159, :cnt 5} {:sz 160, :cnt 1} {:sz 162, :cnt 1} {:sz 163, :cnt 1} {:sz 186, :cnt 1} {:sz 199, :cnt 1} {:sz 200, :cnt 1} {:sz 205, :cnt 1} {:sz 211, :cnt 1} {:sz 212, :cnt 1} {:sz 213, :cnt 1} {:sz 216, :cnt 1} {:sz 222, :cnt 1} {:sz 223, :cnt 1}]


(io/with-out-writer
  "/store/data/PanClus/Entropy/ref77+PG30-SP-clusters-centers-jsd4-cksz-50.clj"
  (prn ref77+PG30-SP-clusters))






(->> ref77-SP-clusters vals
     (#(->> % (map :members) (filter #(> (count %) 1))))
     (group-by #(-> % :jsd :m (roundit :places 2)))
     (map (fn[[k m]] [k (count m)]))
     (sort-by first)
     (mapv (fn[[jsd ct]] {:jsd jsd :cnt ct})))
