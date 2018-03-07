(ns panclus.clustering
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
   [aerial.bio.utils.aligners :as aln]
   [aerial.bio.utils.seqxform :as sqx]

   [panclus.params :as pams]
   ])




(defn roundit [r & {:keys [places] :or {places 4}}]
  (let [n (Math/pow 10.0 places)]
    (-> r (* n) Math/round (/ n))))




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

(defn member->sq
  "Return the sequence for member id 'member' which is
  entry,locus-tag. 'aa' indicates whether or not the sq should be the
  amino acid translation."
  [member & {:keys [aa] :or {aa true}}]
  (let [xform (if aa ntseq2aaseq identity)
        x (member->entry member)
        genomebase (pams/get-params [:biodb-info :genomes :base])
        basedir (if (= (str/substring x 0 3) "TVO")
                  (fs/join genomebase (pams/get-params
                                       [:biodb-info :genomes :tvoseq02]))
                  (fs/join genomebase (pams/get-params
                                       [:biodb-info :genomes :refseq77])))]
    (->> (bufiles/gen-name-seq x :basedir basedir) second xform)))

(defn member->dist
  "Return the kmer distribution for member id 'member' which is
  entry,locus-tag. 'aa' indicates whether or not the sq used should be
  the ammino acide translation. It also indicates the word size off of
  var 'ows'"
  [member & {:keys [aa] :or {aa true}}]
  (let [sq (member->sq member :aa aa)
        ws (ows (if aa :aa :nt))]
    (p/probs ws sq)))

(defn member->clus-input
  [member {:keys [aa] :or {aa true}}]
  (let [sq (member->sq member :aa aa)
        cnt (count sq)
        dist (p/probs (ows (if aa :aa :nt)) sq)]
      [member cnt dist]))


(defn get-dist-cnt [x & {:keys [aa] :or {aa true}}]
  (if (map? x)
    x
    (let [xsq (member->sq x)
          xcnt (count xsq)
          xdist (p/probs (ows (if aa :aa :nt)) xsq)]
      [xdist xcnt])))

(defn get-dist [x & {:keys [aa] :or {aa true}}]
  (if (map? x)
    x
    (let [[dist cnt] (get-dist-cnt x :aa aa)]
      dist)))


(defn id-line->member
  [idl]
  (let [[ent & x] (->> idl (str/split #","))
        lt (->> x (coll/partitionv-all 2) (into {})
                (#(% "locus_tag")))]
    (str ent "," lt)))

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
  (if (string? strain-fna)
    (let [name-seq-pairs (bufiles/read-seqs strain-fna :info :both)]
      (gen-sqs-dists
       (sort-by second (fn[l r] (< (count l) (count r))) name-seq-pairs)
       :aa aa :ows ows))
    (vfold member->clus-input
           strain-fna)))

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


(comment

  (def ref-dists-all (gen-strain-group-dists reffnas :ows 6))
  (def center-dists (map (fn[tuple]
                           {:center tuple :members #{} :name (first tuple)})
                         (first ref-dists-all)))
  (def ref-dists (rest ref-dists-all))
  (def refdists-map (gen-nm-dists-map ref-dists))
  (def ref-len-dists-map (gen-len-dists-map ref-dists))

  (->> ref-dists-all
     (concat (gen-strain-group-dists (pgfnas :pg30) :ows (ows :aa)))
     (apply concat) count)

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

(defn new-clu
  "Return a new cluster start node from initial center element defined
  by [nm cnt Q]. Where nm is a 'entry,locus-tag' id, cnt is the length
  of the sequence associated with nm, and Q is the kmer distribution
  of the sequence.

  Returns a map with entries [k m], with k a PGSP locus tag formed
  using current value of PGSP-index and m is the map:

  {:center [k cnt Q], :dists [Q] :jsds [n 0.0]
  :cnt {:sm cnt, :d 1, :m cnt}
  :members #{n} :name k}"

  [nm cnt Q]
  (let [pgnm (format "PGSP_%05d" (swap! PGSP-index inc))]
    {:center [pgnm cnt Q]
     :dists [Q]
     :jsds [[nm 0.0]]
     :cnt {:sm cnt, :d 1, :m cnt}
     :members #{nm}
     :name pgnm}))

(defn add-member
  "Add a new element 'member' to the current members of the cluster clu
  identified by the center element 'center', which is a vector of the
  form [pgnm jsd cnt Q], where pgnm is a PGSP locus tag. 'clus' is a
  clustering (map of pg locus tags to clusters), where [pgnm clu] is an
  entry."
  [clus center member]
  (let [[n1 pcnt P] member
        [nm jsd cnt Q] center
        curclus (get clus nm)
        members (:members curclus)
        cur-dists (:dists curclus)
        jsds (conj (curclus :jsds) [n1 jsd])
        new-cnt (-> curclus :cnt (running-mean pcnt))]
    (-> curclus
        (assoc-in [:dists] (conj cur-dists P))
        (assoc-in  [:members] (conj members n1))
        (assoc-in [:jsds] jsds)
        (assoc-in [:cnt] new-cnt))))

(defn best-center
  "Compute the cluster whose center best matches the distribution
  given by P. 'ds' is a collection of cluster centers, each of
  standard form [pgnm cnt Q], where pgnm is the PG locus tag for
  cluster, cnt is the mean of current member lengths, and Q is the
  current centroid distribution. Returns a vector [pgnm jsd cnt Q],
  where jsd is JSD score of P to Q. If ds is empty, returns nil."
  [ds P]
  (->> ds
       (vfold (fn[[n cnt Q]]
                [n (roundit (it/jensen-shannon P Q)) cnt Q]))
       (sort-by second)
       first))

(defn candidate-centers
  "Compute the set of centers to be used in comparisons with a
  sequence distribtion. These are biased toward the %-max/%-max+
  length differences to the candidate as given by pcnt. len-distmap is
  a map associating center lengths to distributions"
  [len-distmap pcnt %-max %-max+]
  (->> (range (-> pcnt (* %-max) Math/round)
              (-> pcnt (* %-max+) Math/round))
       (mapcat len-distmap)))


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


(defn clustering->singlet-dists
  "Obtain and return the collection of singlets from 'clustering' in
  the form of std input distribution elements: [name, cnt dist], where
  name is the entry,locus_tag id for sequence, count is sequence
  length, and dist is the associated sequence distribution. All of
  those elements are part of the cluster element for a singlet in
  'clustering'"
  [clustering]
  (->> clustering vals
       (filter #(-> % :members count (= 1)))
       (mapv (fn[m]
               (let [[nm cnt dist] (m :center)
                     name   (-> m :members first)]
                 [name cnt dist])))))

(defn clustering-minus-singlets
  "Obtain the collection of clusters in 'clustering' which are not
  singlets - clusters have > 1 member. Return this collection, which
  can be used as start-centers for a new clustering."
  [clustering]
  (->> clustering
       (filter (fn[[n m]]
                 (-> m :members count (> 1))))
       (into {})))


(defn finish-centers
  ""
  [start-clus %-max %-max+]
  (let [clus (make-hybrids (clustering-minus-singlets start-clus))
        len-distmap (gen-len-dists-map (->> clus vals (mapv #(% :center))))
        singlets (apply dissoc start-clus (keys clus))]
    (reduce (fn[M [nm clu]]
              (let [P (-> clu :center last)
                    pcnt (-> clu :center second)
                    pt [(-> clu :members first) pcnt P]
                    ds (candidate-centers len-distmap pcnt %-max %-max+)
                    best (best-center ds P)]
                (if (and best (<= (second best) 0.2))
                  (assoc M (first best) (add-member clus best pt))
                  (assoc M (clu :name) clu))))
            clus singlets)))

(defn init-centers
  ""
  [dists & {:keys [%-max clus] :or {%-max 0.8}}]
  (let [ndists (gen-nm-dists-map dists)
        %-max+ (->> %-max (- 1.0) (+ 1.0) double)]
    (loop [ndists ndists
           clus (or clus {})]
      #_(println (format "cnt ndists: %s; cnt clus %s"
                       (count ndists) (count clus)))
      (if (empty? ndists)
        (-> clus (finish-centers %-max %-max+) make-hybrids)
        (let [[nm pcnt P :as pt] (-> ndists first second)
              len-distmap (gen-len-dists-map
                           (->> clus vals (mapv #(% :center))))
              ds (->> (range (-> pcnt (* %-max) Math/round)
                             (-> pcnt (* %-max+) Math/round))
                      (mapcat len-distmap))
              best (best-center ds P)
              clu (if (and best (<= (second best) 0.2))
                    (add-member clus best pt)
                    (new-clu nm pcnt P))]
          (recur (dissoc ndists nm)
                 (assoc clus (clu :name) clu)))))))

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
      (init-centers coll))))


(defn cluster-strains
  ""
  [dists start-clusters
   & {:keys [%-max jsdctpt pgsp-start]
      :or {%-max 0.80 jsdctpt 0.40 pgsp-start @PGSP-index}}]

  (let [%-max+ (->> %-max (- 1.0) (+ 1.0) double)
        dists (->> dists (apply concat))
        clus (init-clusters start-clusters pgsp-start)
        center-dists (mapv #(% :center) (vals clus))
        len-distmap (gen-len-dists-map center-dists)]
    (loop [dists dists
           clus clus]
      (if (empty? dists)
        clus
        (let [[n1 pcnt P :as pt] (first dists)
              ;;_ (prn :n1 n1 :pcnt pcnt :means mean-std mean+std)

              ;; Adjust clus distribution set to use, for current
              ;; dists point, by length - Only for optimization - not
              ;; actually needed for RE
              ds (->> (range (-> pcnt (* %-max) Math/round)
                             (-> pcnt (* %-max+) Math/round))
                      (mapcat len-distmap))
              ;; Gen scores of current dists point to centers in 'ds'
              ;; set. Pick center with lowest score.
              best-center (best-center ds P)
              ;; If best-center score <= jsdctpt, add current point to
              ;; the cluster for this center
              clu (when (and best-center (<= (second best-center) jsdctpt))
                    (add-member clus best-center pt))
              new-clus (if clu (assoc clus (clu :name) clu) clus)]
          #_(println :new-distmap (count new-distmap))
          (recur (rest dists)
                 new-clus))))))


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
  clusters is a map with keys [center members dists :cnt :jsds :name]
  and each member of start-dists is standard tuple [n cnt dist].
  Returns collection of tuples in start-dists not in any cluster."
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
   & {:keys [%-max jsdctpt diffcut pgsp-start]
      :or {%-max 0.8 jsdctpt 0.40 diffcut 0}}]
  (let [all-dists (apply concat dists)]
    (loop [clusters (cluster-strains
                     dists center-dists
                     :%-max %-max :jsdctpt jsdctpt
                     :pgsp-start pgsp-start)
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
  [dists & {:keys [%-max jsdctpt]
            :or {%-max 0.8 jsdctpt 0.40}}]
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
                     :%-max %-max :jsdctpt jsdctpt
                     :pgsp-start @PGSP-index)]
      (merge-clusters
       new-clus
       (cluster-leftovers (leftovers-from-clustering new-clus dists)
                          :%-max %-max :jsdctpt jsdctpt)))))


(defn final-pass
  "At the end of a strain clustering, this is a final pass over the
  current results, as embodied in 'clustering', to attempt singlet
  coalescing and/or inclusion in current multi-member clusters. This
  final pass uses a %-max of 0.8 to check final candidates against
  centers with up to a 40% length difference and jsdctpt of 0.4,
  irrespective of the original values given for the strain
  clustering."
  [clustering]
  (let [singlet-dists (clustering->singlet-dists clustering)
        start-centers (clustering-minus-singlets clustering)
        pgsp-num (->> clustering keys sort last
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
   & {:keys [%-max jsdctpt diffcut chunk-size pgsp-start]
      :or {%-max 0.8 jsdctpt 0.40 diffcut 0 chunk-size 10}}]
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
                              :%-max %-max :jsdctpt jsdctpt
                              :pgsp-start pgsp-num :diffcut diffcut)

              leftover-clus (cluster-leftovers
                             (leftovers-from-clustering
                              chunk-clusters (apply concat chunk-dists))
                             :%-max %-max :jsdctpt jsdctpt)]
          (println "At recur point")
          (recur @PGSP-index
                 (drop chunk-size strain-fnas)
                 (merge-clusters chunk-clusters leftover-clus)))))))


(defn pair-jsd [x y & {:keys [aa] :or {aa true}}]
  (let [xdist (get-dist x :aa aa)
        ydist (get-dist y :aa aa)]
    (it/jensen-shannon xdist ydist)))




(defn center-dist->pair-map
  [center-dist jsdct & {:keys [dir] :or {dir :below}}]
  (let [func (if (= dir :below) coll/take-until coll/drop-until)]
    (->> center-dist
         (sort-by first)
         (func #(-> % first (> jsdct)))
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
                 {}))))

(defn next-set
  "Transitive merge set computation. 'curset' is the current total set
  of related clusters to merge, 'pair-map' is the current map of
  leftover possible clusters to add in, 'mems' is the current
  transitive elements to be checked. Returns [curset curmap], where
  curset here is the new set of mergeable clusters (may be the same as
  input curset), curmap is the current version of pair-map with any
  newly added cluster keys removed."
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

(defn next-id
  "Generate a new clu-id for merged cluster set"
  []
  (format "PGSP_%05d" (swap! PGSP-index inc)))

(defn tran-reduce
  "Driver for generating the set of mergable clusters from an
  originating 'pair-map'. 'pair-map' is the map of all pairwise
  clusters satisfying a maximal center-center JSD cutoff. 'pgsp-start'
  is an integer to start generating new cluster ids for the mergable
  cluster sets. Generally it will be the value of PGSP-index after the
  clustering being used to generate merge sets."
  [pair-map pgsp-start]
  (swap! PGSP-index (constantly pgsp-start))
  (reduce (fn[[M newmap] pgk]
            (if (not (newmap pgk))
              [M newmap]
              (let [[cluset newmap] (next-set #{} newmap [pgk])
                    clu-id (next-id)]
                [(assoc M clu-id cluset) newmap])))
          [{} pair-map] (keys pair-map)))


(defn merged-set->cluster
  [merged-set clustering & {:keys [ows] :or {ows 6}}]
  (let [[pgnm pgids] merged-set
        mems (->> pgids
                  (mapv #(->> (clustering %) :members))
                  (apply set/union))

        oldoutjsds (->> pgids
                     (mapv #(->> (clustering %) :outjsds))
                     (apply coll/concatv))
        oldjsds (->> pgids
                     (mapv #(->> (clustering %) :jsds))
                     (mapv #(sort-by second %)))
        bestjsds (mapv first oldjsds)
        bestdists (mapv #(-> % first member->dist) bestjsds)
        centroid (it/hybrid-dictionary ows bestdists)

        mdc-trips (mapv #(let [ent (member->entry %)
                               [dist cnt] (get-dist-cnt %)]
                           [% dist cnt])
                        mems)
        jsds (->> mdc-trips
                  (mapv (fn[[mem P _]]
                          [mem (roundit (it/jensen-shannon P centroid))]))
                  (sort-by second))
        outjsds (coll/drop-until #(-> % second (>= 0.31)) jsds)
        jsds (take-while #(-> % second (< 0.31)) jsds)
        gdmems (->> jsds (mapv first) (into #{}))
        mdc-trips (filter #(-> % first gdmems) mdc-trips)

        cnt (->> mdc-trips (mapv last) p/mean)
        Q (->> mdc-trips (mapv second) (it/hybrid-dictionary ows))

        jsds (mapv (fn[[mem P _]] [mem (roundit (it/jensen-shannon P Q))])
                   mdc-trips)]

    {:center [pgnm cnt Q]
     :dists [Q]
     :jsds jsds
     :outjsds (coll/concatv outjsds oldoutjsds)
     :cnt {:sm cnt, :d 1, :m cnt}
     :members gdmems
     :name pgnm}))

(defn merged-sets->clusters
  [merged-sets clustering]
  (->> merged-sets
       (vfold (fn[ms]
                (let [clu (merged-set->cluster ms clustering)]
                  [(clu :name) clu])))
       (into {})))


(def DBG (atom {}))

(defn cur->next-merged-clusters
  [cur-clusters pgsp-start pair-map]
  (let [next-merged-sets (-> pair-map
                             (tran-reduce pgsp-start) first)
        next-merged-clusters (merged-sets->clusters
                              next-merged-sets cur-clusters)]
    (swap! DBG (fn[m] (assoc m :next-merged-sets next-merged-sets
                              :next-merged-clusters next-merged-clusters)))
    next-merged-clusters))

(defn cur->next-clustering
  [cur-clusters pair-map next-merged-clusters]
  (let [cur-clus-minus (->> pair-map keys (apply dissoc cur-clusters))
        cur-merged-clus (merge cur-clus-minus next-merged-clusters)]
    (swap! DBG
           (fn[m] (assoc m :current-clustering-minus cur-clus-minus
                          :current-merged-clustering cur-merged-clus)))
    cur-merged-clus))


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

(defn center-distances
  [cur-clusters next-merged-clusters cur-center-distances jsdct]
  (let [center-jsd-dist cur-center-distances
        new-keys (-> next-merged-clusters keys set)
        cur-keys (-> cur-clusters keys set)]
    (->> center-jsd-dist (sort-by first)
         (coll/drop-until #(-> % first (> jsdct)))
         (mapcat (fn[[scr members]] (mapcat identity members)))
         (into #{})
         (set/difference cur-keys)
         (mapv #(vector % (cur-clusters %)))
         (into {})
         vals
         compare-clusters)))


(defn cur-members-from-distances
  [cur-center-distances pair-map jsdct]
  (let [above-jsdct (->> cur-center-distances
                         (sort-by first)
                         (coll/drop-until #(-> % first (> jsdct)))
                         (mapcat (fn[[scr members]] (mapcat identity members)))
                         (into #{})
                         (set/intersection (-> pair-map keys set))
                         (mapcat #(conj (pair-map %) %))
                         (into #{}))]))

(comment


(def tigr4-centers
  (init-clusters (-> reffnas first (gen-strain-dists :ows (ows :aa))) 0))

(io/with-out-writer "/store/data/PanClus/Entropy/tigr4-centers.clj"
  (prn tigr4-centers))

(def tigr4-centers
  (-> (pams/get-params :panclus-base)
      (fs/join "Entropy/tigr4-centers.clj")
      slurp read-string))



(swap! PGSP-index
       (fn[& args]
         (->> tigr4-centers
              keys sort last (str/split #"_") second Integer.)))

(def ref77-SP-clusters
  (future (run-strain-clustering
           (rest reffnas)
           tigr4-centers
           :pgsp-start @PGSP-index :jsdctpt 0.4 :%-max 0.80
           :chunk-size 2 :diffcut 0)))

(def ref77-SP-clusters @ref77-SP-clusters)

(io/with-out-writer "/store/data/PanClus/Entropy/ref77-SP-clusters.clj"
  (prn ref77-SP-clusters))

(def ref77-SP-clusters
  (-> (pams/get-params :panclus-base)
      (fs/join "Entropy/ref77-SP-clusters.clj")
      slurp read-string))

(swap! PGSP-index
       (fn[& args]
         (->> ref77-SP-clusters
              keys sort last (str/split #"_") second Integer.)))



(def ref77+PG30-SP-clusters
  (future (run-strain-clustering
           (pgfnas :pg30)
           ref77-SP-clusters
           :pgsp-start @PGSP-index :jsdctpt 0.4 :%-max 0.80
           :chunk-size 2 :diffcut 0)))

(def ref77+PG30-SP-clusters @ref77+PG30-SP-clusters)

(io/with-out-writer "/store/data/PanClus/Entropy/ref77+PG30-SP-clusters.clj"
  (prn ref77+PG30-SP-clusters))

(def ref77+PG30-SP-clusters
  (-> (pams/get-params :panclus-base)
      (fs/join "Entropy/ref77+PG30-SP-clusters.clj")
      slurp read-string))

(swap! PGSP-index
       (fn[& args]
         (->> ref77+PG30-SP-clusters
              keys sort last (str/split #"_") second Integer.)))




(def ref77pg350-fut
  (future (run-strain-clustering
           (pgfnas :pg320)
           ref77+PG30-SP-clusters
           :pgsp-start @PGSP-index :jsdctpt 0.4 :%-max 0.95
           :chunk-size 2 :diffcut 0)))

(io/with-out-writer "/store/data/PanClus/Entropy/ref77+PG350-SP-clusters.clj"
  (prn @ref77pg350-fut))

(def ref77pg350-fut
  (future
    (-> (pams/get-params :panclus-base)
        (fs/join "Entropy/ref77+PG350-SP-clusters.clj")
        slurp read-string)))

(swap! PGSP-index
       (fn[& args]
         (->> @ref77pg350-fut
              keys sort last (str/split #"_") second Integer.)))




;;; full ref77pg350 merged clusters
(def merged-clusters
  (merged-sets->clusters merged-sets))

(io/with-out-writer
  "/store/data/PanClus/Entropy/full-r77pg350-merged-clusters.clj"
  (prn merged-clusters))

(def merged-clusters
  (-> (pams/get-params :panclus-base)
      (fs/join "Entropy/full-r77pg350-merged-clusters.clj")
      slurp read-string))


;;; ref77pg350 minus all input clusters to merge plus new merged clusters
(def ref77pg350-minus
  (->> testmap keys (apply dissoc @ref77pg350-fut) count))

(def ref77-pg350-with-merged-clusters
  (merge ref77pg350-minus merged-clusters))

(io/with-out-writer
  "/store/data/PanClus/Entropy/ref77-pg350-with-merged-clusters.clj"
  (prn ref77-pg350-with-merged-clusters))

(def ref77-pg350-with-merged-clusters
  (-> (pams/get-params :panclus-base)
      (fs/join "Entropy/ref77-pg350-with-merged-clusters.clj")
      slurp read-string))

(->> ref77-pg350-with-merged-clusters
     keys sort last (str/split #"_") second Integer.)





)
