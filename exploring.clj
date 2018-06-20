
;;; Pair down # seqs to compare...


(def len-map (len-map-data))
  #_(-> (pams/get-params :panclus-base)
      (fs/join "Data/all-strain-gene-aa-lens.clj")
      slurp read-string)


(def _19F-22F-mems
  (->> _19F-22F-fnas
       (vfold #(bufiles/read-seqs % :info :name))
       (apply concat)
       (map (fn[n]
              (let [[ent & x] (->> n (str/split #","))
                    lt (->> x (coll/partitionv-all 2) (into {})
                            (#(% "locus_tag")))]
                (str ent "," lt))))))

(def _19F-22F-ident-map
  (let [len-map (merge len-map jsa)
        m (->> _19F-22F-mems
               (mapv #(vector % (len-map %)))
               (group-by second)
               (map (fn[[k mems]] [k (mapv first mems)])))]
    (->>
     (mapcat
      (fn[[cnt mems]]
        (loop [mems (->> mems (vfold #(vector % (member->sq %))) (into {}))
               grps []]
          (println (count mems) (map #(-> % second count) grps))
          (if (= (count mems) 0)
            grps
            (let [xnm (-> mems first first)
                  x (-> mems first second
                        #_(#(do (println :mem (first mems) :cnt (count %)) %))
                        (str/substring 2 -2))
                  ys (rest mems)
                  res (->> ys
                           (vfold (fn[[nm y]] [nm (->> y (str/substring? x))]))
                           (group-by second)
                           (map (fn[[k mems]] [k (mapv first mems)]))
                           (into {}) (#(% true)))
                  cnt (count res)]
              (if (> cnt 0)
                (recur (apply dissoc mems (conj res xnm))
                       (conj grps [xnm res]))
                (recur (dissoc mems xnm)
                       (conj grps [xnm [xnm]])))))))
      m)
     (into {}))))


(def r77mems
  (->> reffnas
       (vfold #(bufiles/read-seqs % :info :name))
       (apply concat)
       (map (fn[n]
              (let [[ent & x] (->> n (str/split #","))
                    lt (->> x (coll/partitionv-all 2) (into {})
                            (#(% "locus_tag")))]
                (str ent "," lt))))))

(->> r77mems (coll/drop-until #(-> % member->strain (not= "NC_003028")))
     count)
(->> r77mems (mapv #(vector % (len-map %))) (into {}) count)

(def r77-ident-map
  (let [m (->> r77mems
               (coll/drop-until #(-> % member->strain (not= "NC_003028")))
               (mapv #(vector % (len-map %)))
               (group-by second)
               (map (fn[[k mems]] [k (mapv first mems)])))]
    (->>
     (mapcat
      (fn[[cnt mems]]
        (loop [mems (->> mems (vfold #(vector % (member->sq %))) (into {}))
               grps []]
          (println (count mems) (map #(-> % second count) grps))
          (if (= (count mems) 0)
            grps
            (let [xnm (-> mems first first)
                  x (-> mems first second (str/substring 2 -2))
                  ys (rest mems)
                  res (->> ys
                           (vfold (fn[[nm y]] [nm (->> y (str/substring? x))]))
                           (group-by second)
                           (map (fn[[k mems]] [k (mapv first mems)]))
                           (into {}) (#(% true)))
                  cnt (count res)]
              (if (> cnt 0)
                (recur (apply dissoc mems (conj res xnm))
                       (conj grps [xnm res]))
                (recur (dissoc mems xnm)
                       (conj grps [xnm [xnm]])))))))
      m)
     (into {}))))

(io/with-out-writer "/store/data/PanClus/Data/r77-seqs-equal-map.clj"
  (prn r77-ident-map))
(def r77-ident-map
  (-> (pams/get-params :panclus-base) (fs/join "Data/r77-seqs-equal-map.clj")
      slurp read-string))

(io/with-out-writer "/store/data/PanClus/Stats/r77-seqs-len-dist.clj"
  (prn (->> r77-ident-map
            (mapv (fn[[k v]] {:len (len-map k) :cnt (inc (count v))})))))

[18067 45570 0.3964669738863287]




(def pg30mems
  (->> (pgfnas :pg30)
       (vfold #(bufiles/read-seqs % :info :name))
       (apply concat)
       (map (fn[n]
              (let [[ent & x] (->> n (str/split #","))
                    lt (->> x (coll/partitionv-all 2) (into {})
                            (#(% "locus_tag")))]
                (str ent "," lt))))))

(->> pg30mems (mapv #(vector % (len-map %))) (into {}) count)

(def pg30-ident-map
  (let [m (->> pg30mems
               (mapv #(vector % (len-map %)))
               (group-by second)
               (map (fn[[k mems]] [k (mapv first mems)])))]
    (->>
     (mapcat
      (fn[[cnt mems]]
        (loop [mems (->> mems (vfold #(vector % (member->sq %))) (into {}))
               grps []]
          (println (count mems) (map #(-> % second count) grps))
          (if (= (count mems) 0)
            grps
            (let [xnm (-> mems first first)
                  x (-> mems first second (str/substring 2 -2))
                  ys (rest mems)
                  res (->> ys
                           (vfold (fn[[nm y]] [nm (->> y (str/substring? x))]))
                           (group-by second)
                           (map (fn[[k mems]] [k (mapv first mems)]))
                           (into {}) (#(% true)))
                  cnt (count res)]
              (if (> cnt 0)
                (recur (apply dissoc mems (conj res xnm))
                       (conj grps [xnm res]))
                (recur (dissoc mems xnm)
                       (conj grps [xnm [xnm]])))))))
      m)
     (into {}))))

(io/with-out-writer "/store/data/PanClus/Data/pg30-seqs-equal-map.clj"
  (prn pg30-ident-map))
(def pg30-ident-map
  (-> (pams/get-params :panclus-base) (fs/join "Data/pg30-seqs-equal-map.clj")
      slurp read-string))

(io/with-out-writer "/store/data/PanClus/Stats/pg30-seqs-len-dist.clj"
  (prn (->> pg30-ident-map
            (mapv (fn[[k v]] {:len (len-map k) :cnt (inc (count v))})))))

[19775 64769 0.3053158146644228]




(def pg320-ident-map
  (let [m (->> (apply dissoc len-map (concat r77mems pg30mems))
               (group-by second)
               (map (fn[[k mems]] [k (mapv first mems)])))]
    (->>
     (mapcat
      (fn[[cnt mems]]
        (loop [mems (->> mems (vfold #(vector % (member->sq %))) (into {}))
               grps []]
          (println (count mems) (map #(-> % second count) grps))
          (if (= (count mems) 0)
            grps
            (let [xnm (-> mems first first)
                  x (-> mems first second (str/substring 2 -2))
                  ys (rest mems)
                  res (->> ys
                           (vfold (fn[[nm y]] [nm (->> y (str/substring? x))]))
                           (group-by second)
                           (map (fn[[k mems]] [k (mapv first mems)]))
                           (into {}) (#(% true)))
                  cnt (count res)]
              (if (> cnt 0)
                (recur (apply dissoc mems (conj res xnm))
                       (conj grps [xnm res]))
                (recur (dissoc mems xnm)
                       (conj grps [xnm [xnm]])))))))
      m)
     (into {}))))

(io/with-out-writer "/store/data/PanClus/Data/pg320-seqs-equal-map.clj"
  (prn pg320-ident-map))
(def pg320-ident-map
  (-> (pams/get-params :panclus-base) (fs/join "Data/pg320-seqs-equal-map.clj")
      slurp read-string))

(io/with-out-writer "/store/data/PanClus/Stats/pg320-seqs-len-dist.clj"
  (prn (->> pg320-ident-map
            (mapv (fn[[k v]] {:len (len-map k) :cnt (inc (count v))})))))

[50680 692167 0.07321932423822575]



(let[memcnt (->> pg320-ident-map vals (map set) (apply set/union) count)
     keycnt (-> pg320-ident-map keys set count)
     ovlapcnt (count (set/intersection
                      (-> pg320-ident-map keys set)
                      (->> pg320-ident-map vals (map set) (apply set/union))))
     totcnt (+ memcnt keycnt (- ovlapcnt))]
  [keycnt totcnt (-> keycnt (/ totcnt) double)])

