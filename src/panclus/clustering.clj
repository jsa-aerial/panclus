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
;;; (separate-strains "refseq77-genes-locus-tags-t4-pa-ab-ec.csv"
;;;                   "RefSeq77/Xspecies")
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
;;; (future  (gen-strain-fastas "PG/PG320" :tvoseq02))
;;; (gen-strain-fastas "RefSeq77/Xspecies" :refseq77)
(defn gen-strain-fastas [strain-dir db]
  (let [base panclus-base
        strainbase (fs/join base strain-dir)
        tags ["gene" "old_locus_tag" "locus_tag" "protein_id"]
        csvs (fs/re-directory-files (fs/join strainbase "CSV") "*.csv")
        dbresolver (aerial.bio.utils.params/get-params [:genomes db])]
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


;;; T4, Pseudomonas, Ecoli blessed, Ecoli old, Abuaumanii
(def species-fnas
  (->> (fs/glob (fs/join panclus-base "RefSeq77/Xspecies/N*.fna"))
       sort vec))

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


(def ^:dynamic ows {:nt 9 :aa 6})

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
  [member & {:keys [aa] :or {aa true}}]
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


(defn pair-jsd [x y & {:keys [aa] :or {aa true}}]
  (let [xdist (get-dist x :aa aa)
        ydist (get-dist y :aa aa)]
    (it/jensen-shannon xdist ydist)))


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
                [n (roundit (it/jensen-shannon P Q) :places 2) cnt Q]))
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
  [clusters & {:keys [wz] :or {wz (ows :aa)}}]
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

(defn calc-max+
  "Calculate the max+ for a sequence length comparison from %-max
  divergence from test sequence. %-max is in open interval (0.0
  .. 1.0). If %-max >= 0.05, %-max+ computed as divergence from test
  seq length as percentage:

  (->> %-max (- 1.0) (+ 1.0) double).

  If %-max < 0.05, switch to multiples of test sequence length:

  (->> %-max (/ 1.0) double)"
  [%-max]
  (if (< %-max 0.05) ; hackish way to switch divergence
    (->> %-max (/ 1.0) double)
    (->> %-max (- 1.0) (+ 1.0) double)))

(defn init-centers
  ""
  [dists & {:keys [%-max clus] :or {%-max 0.4}}]
  (let [ndists (gen-nm-dists-map dists)
        %-max+ (calc-max+ %-max)]
    (loop [ndists ndists
           clus (or clus {})]
      #_(println (format "cnt ndists: %s; cnt clus %s"
                       (count ndists) (count clus)))
      (if (empty? ndists)
        (-> clus (finish-centers %-max %-max+) make-hybrids)
        (let [[nm pcnt P :as pt] (-> ndists first second)
              len-distmap (gen-len-dists-map
                           (->> clus vals (mapv #(% :center))))
              ds (candidate-centers len-distmap pcnt %-max %-max+)
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

  (let [%-max+ (calc-max+ %-max)
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
              ds (candidate-centers len-distmap pcnt %-max %-max+)

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
    (swap! PGSP-index (constantly pgsp-start))
    (loop [dists dists
           centers center-dists
           leftcnt (count all-dists)]
      #_(prn "         (first clusters)" (first clusters))
      (let [clusters (cluster-strains
                      dists centers
                      :%-max %-max :jsdctpt jsdctpt)
            leftovers (leftovers-from-clustering clusters all-dists)
            cur-leftcnt (count leftovers)
            newclus (make-hybrids clusters)]
        (println :diff (- leftcnt cur-leftcnt))
        (if (<= (- leftcnt cur-leftcnt) diffcut)
          newclus
          (recur [leftovers]
                 newclus
                 cur-leftcnt))))))


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


(defn write-cluster [clus name]
  (io/with-out-writer
    (-> (pams/get-params :panclus-base) (fs/join "Data" name))
    (prn clus))
  clus)

(defn final-pass
  "At the end of a strain clustering, this is a final pass over the
  current results, as embodied in 'clustering', to attempt singlet
  coalescing and/or inclusion in current multi-member clusters. This
  final pass uses a %-max of 0.8 to check final candidates against
  centers with up to a 40% length difference and jsdctpt of 0.4,
  irrespective of the original values given for the strain
  clustering."
  [clustering jsdctpt %-max filename]
  (let [singlet-dists (clustering->singlet-dists clustering)
        start-centers (clustering-minus-singlets clustering)
        pgsp-num (->> clustering keys sort last
                      (str/split #"_") second Integer.)
        newclus (run-clustering-cycle
                 [singlet-dists] start-centers
                 :%-max %-max :jsdctpt jsdctpt
                 :pgsp-start pgsp-num)]
    (write-cluster
     (merge-clusters
      newclus
      (cluster-leftovers
       (leftovers-from-clustering newclus singlet-dists) :%-max %-max))
     filename)))


(defn run-strain-clustering
  ""
  [strain-fnas centers-start
   & {:keys [%-max final-%-max jsdctpt final-jsdctpt
             diffcut chunk-size pgsp-start]
      :or {%-max 0.8 final-%-max 0.1 jsdctpt 0.40 final-jsdctpt 0.4
           diffcut 0 chunk-size 10}}]
  (println :let)
  (let [alpha :aa
        ows (ows alpha)]
    (println :at-loop)
    (loop [pgsp-num pgsp-start
           strain-fnas strain-fnas
           names (cycle ["snapshot-clus-1.clj"
                         "snapshot-clus-2.clj"])
           clusters centers-start]
      (if (empty? strain-fnas)
        (final-pass
         clusters final-jsdctpt final-%-max "final-snapshot-clus.clj")

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
                             :%-max %-max :jsdctpt jsdctpt)
              new-clus (merge-clusters chunk-clusters leftover-clus)]

          ;; Snapshot current results - fail over and other fuckups
          (write-cluster new-clus (first names))

          (println "At recur point")
          (recur @PGSP-index
                 (drop chunk-size strain-fnas)
                 (rest names)
                 new-clus))))))


(defn min2full-clustering
  [min-clustering mem-ident-map]
  (->> min-clustering
       (map (fn[[nm m]]
              (let [mems (m :members)
                    fullmems (->> mems
                                  (mapcat #(mem-ident-map %))
                                  (into #{})
                                  (set/union mems))
                    jsds (m :jsds)
                    fulljsds (->> jsds
                                  (mapcat (fn[[nm jsd]]
                                            (->> nm
                                                 mem-ident-map
                                                 (#(conj % nm)) set
                                                 (map #(vector % jsd)))))
                                  vec)
                    newm (assoc m :members fullmems
                                  :jsds fulljsds)]
                [nm newm])))
       (into {})))



(defn next-point-dist [clusters mode]
  (if (= mode :center-center)
    (-> clusters first last)
    (assert :NYI "only :center-center supported at this point")))

(defn get-comp-pairs [dists Pnm P]
  (->> dists
       (vfold (fn[[nm Q]]
                [(roundit (it/jensen-shannon P Q) :places 2) [Pnm nm]]))))

(defn new-comp-dist [comp-dist comp-pairs]
  (reduce (fn[CD [scr pair]]
            (assoc CD scr (conj (get CD scr []) pair)))
          comp-dist comp-pairs))

(defn compare-clusters [clusters & {:keys [mode] :or {mode :center-center}}]
  (loop [comp-dist {}
         clusters (mapv #(vector (% :name) (-> % :center last)) clusters)]
    (if (empty? (rest clusters))
      comp-dist
      (let [P (next-point-dist clusters mode)
            Pnm (-> clusters first first)
            dists (rest clusters)]
        (println (format "(count clusters): %s" (count clusters)))
        (recur (new-comp-dist comp-dist (get-comp-pairs dists Pnm P))
               (rest clusters))))))



(def len-map nil)
(defn len-map-data []
  (-> (pams/get-params :panclus-base)
      (fs/join "Data/all-strain-gene-aa-lens.clj")
      slurp read-string))

(defn clus-sqlen-max-min
  [cluid clustering]
  (->> cluid clustering :members
       ((fn[mems] (if (resolve 'len-map)
                   (mapv len-map mems)
                   (vfold #(-> % member->sq count) mems))))
       (#(let [cnt (count %)
               mn (apply min %)
               mx (apply max %)
               ratio (-> mn (/ mx) double (roundit :places 2))]
           {:id cluid :cnt cnt :min mn :max mx :ratio ratio}))))


(defn next-pass-data-sets
  [clustering center-jsd-dist jsdctpt]
  (let [cluids (->>  center-jsd-dist
                     (sort-by first)
                     (coll/take-until #(-> % first (> jsdctpt)))
                     (mapcat (fn[[_ pairs]]
                               (mapcat (fn[pair]
                                         (let [clus (mapv #(clustering %) pair)
                                               nms (mapv #(% :name) clus)]
                                           nms))
                                       pairs)))
                     (into #{}))
        mems (mapcat #(-> % clustering :members) cluids)
        clustering-minus (apply dissoc clustering cluids)]
    [mems clustering-minus cluids]))




(comment


(def tigr4-centers
  (init-clusters (-> reffnas first (gen-strain-dists :ows (ows :aa))) 0))

(io/with-out-writer "/store/data/PanClus/Data/tigr4-centers.clj"
  (prn tigr4-centers))

(def tigr4-centers
  (-> (pams/get-params :panclus-base)
      (fs/join "Data/tigr4-centers.clj")
      slurp read-string))

(swap! PGSP-index
       (fn[& args]
         (->> tigr4-centers
              keys sort last (str/split #"_") second Integer.)))


(def species-clusters
  (future (binding [ows {:nt 9 :aa 5}]
            (run-strain-clustering
             (->> species-fnas rest (take 1))
             (init-clusters (-> reffnas first
                                (gen-strain-dists :ows (ows :aa))) 0)
             :pgsp-start @PGSP-index :jsdctpt 0.37 :%-max 0.2
             :chunk-size 2 :diffcut 0))))


;;; Xref for Stephen
(def _19F-22F-fnas
  ["/store/data/PanClus/RefSeq77/NC_012469.fna"
   "/store/data/PanClus/TVO_1902599.fna"])
(def t4-19f-22f
  (future (run-strain-clustering
           #_(rest reffnas)
           (->> _19F-22F-ident-map keys (partition-all 2000))
           tigr4-centers
           :pgsp-start @PGSP-index :jsdctpt 0.25 :%-max 0.1
           :final-jsdctpt 0.3 :final-%-max 0.04
           :chunk-size 2 :diffcut 0)))



(def ref77-SP-clusters
  (future (run-strain-clustering
           #_(rest reffnas)
           (->> r77-ident-map keys (partition-all 4000))
           tigr4-centers
           :pgsp-start @PGSP-index :jsdctpt 0.3 :%-max 0.3
           :chunk-size 2 :diffcut 0)))

(def ref77-SP-clusters @ref77-SP-clusters)


(io/with-out-writer "/store/data/PanClus/Data/ref77-SP-clusters.clj"
  (prn ref77-SP-clusters))
(io/with-out-writer "/store/data/PanClus/Data/ref77-SP-clusters-done.clj"
  (prn ref77-full-clusters))


(def ref77-SP-clusters
  (-> (pams/get-params :panclus-base)
      (fs/join "Data/ref77-SP-clusters.clj")
      slurp read-string))
(def ref77-SP-clusters
  (-> (pams/get-params :panclus-base)
      (fs/join "Data/EqRes/ref77-SP-clusters-new2.clj")
      slurp read-string))

(swap! PGSP-index
       (fn[& args]
         (->> ref77-SP-clusters
              keys sort last (str/split #"_") second Integer.)))




(def ref77+PG30-SP-clusters
  (future (run-strain-clustering
           #_(pgfnas :pg30)
           (->> pg30-ident-map keys (partition-all 4000))
           ref77-SP-clusters
           :pgsp-start @PGSP-index :jsdctpt 0.3 :%-max 0.3
           :chunk-size 2 :diffcut 0)))

(def ref77+PG30-SP-clusters @ref77+PG30-SP-clusters)

(io/with-out-writer "/store/data/PanClus/Data/EqRes/ref77+PG30-SP-clusters.clj"
  (prn ref77+PG30-SP-clusters))
(io/with-out-writer "/store/data/PanClus/Data/ref77+PG30-SP-clusters-done.clj"
  (prn ref77+pg30-full-clusters))

(def ref77+PG30-SP-clusters
  (-> (pams/get-params :panclus-base)
      (fs/join "Data/ref77+PG30-SP-clusters.clj")
      slurp read-string))

(swap! PGSP-index
       (fn[& args]
         (->> ref77+PG30-SP-clusters
              keys sort last (str/split #"_") second Integer.)))




(def ref77pg350-SP-clusters
  (future (run-strain-clustering
           #_(pgfnas :pg320)
           (->> pg320-ident-map keys (partition-all 4000))
           ref77+PG30-SP-clusters
           :pgsp-start @PGSP-index :jsdctpt 0.3 :%-max 0.3
           :chunk-size 2 :diffcut 0)))

(def ref77pg350-SP-clusters @ref77pg350-SP-clusters)

(io/with-out-writer "/store/data/PanClus/Data/ref77+PG350-SP-clusters.clj"
  (prn ref77pg350-SP-clusters))

(def ref77pg350-SP-clusters
  (-> (pams/get-params :panclus-base)
      (fs/join "Data/ref77+PG350-SP-clusters-1.clj")
      slurp read-string))

(swap! PGSP-index
       (fn[& args]
         (->> ref77pg350-SP-clusters
              keys sort last (str/split #"_") second Integer.)))





)
