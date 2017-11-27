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


