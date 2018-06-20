(ns panclus.core
  "Pangenome clustering server...."

  {:author "Jon Anthony"}

  (:gen-class)

  (:require [clojure.string :as cstr]
            [clojure.tools.reader.edn :as edn]

            ;; Command line arg processing
            [clojure.tools.cli :refer [parse-opts]]

            ;; Tunneling Cider nREPL support
            [clojure.tools.nrepl.server :as nrs]
            #_[cider.nrepl :refer [cider-middleware]] ; BROKEN SIDE EFFECTS
            [refactor-nrepl.middleware :refer [wrap-refactor]]

            [aerial.fs :as fs]
            [aerial.utils.misc :refer [getenv]]
            [aerial.utils.coll :as coll]

            [aerial.bio.utils.params :as bpams]

            [panclus.params :as pams]
            [panclus.clustering] ; just to force loading
            ))

(def cli-options

  [["-r" "--repl-port PORT" "Port for nrepl server"
    :default nil
    :parse-fn #(Integer/parseInt %)
    :validate [#(< 4000 % 5000) "Must be a number between 4001 and 4999"]]
   ;; A boolean option defaulting to nil
   ["-h" "--help" "Print this help"
    :default true]])


;; Incorporate correct cider middleware for tunneling nREPL support
#_(apply nrs/default-handler (map resolve cider-middleware))


(defn- set-configuration []
  (let [m (->> "config.clj" (fs/join (fs/pwd))
               slurp edn/read-string)]
    (swap! pams/params (fn[_] m))
    (bpams/set-configuration (pams/get-params :biodb-info))
    m))


(defn nrepl-handler-hack []
  (require 'cider.nrepl)
  (let [cm (var-get (ns-resolve 'cider.nrepl 'cider-middleware))
        ;;cm (filter #(not= % 'cider.nrepl/wrap-pprint-fn) cm)
        resolve-or-fail (fn[sym] (println :sym sym)
                          (or (ns-resolve 'cider.nrepl sym)
                              (throw (IllegalArgumentException.
                                      (format "Cannot resolve %s" sym)))))]
    (clojure.pprint/pprint cm)
    (apply nrs/default-handler
           (concat (mapv resolve-or-fail cm)
                   [#'wrap-refactor]))))
;;;
#_(apply nrs/default-handler
       (concat (map resolve cider-middleware)
               [#'wrap-refactor]))

(defn -main
  "Daemon starter ...."
  [& args]

  (let [opts (parse-opts args cli-options)
        options (opts :options)
        arguments (opts :arguments)
        summary (opts :summary)
        errors (opts :errors)
        rpl-port (options :repl-port)
        nrepl-handler (nrepl-handler-hack)]
    (if errors
      (do (println "Error(s): " errors)
          (System/exit 1))
      (cond

        rpl-port
        (do
          (println :rpl-port rpl-port)
          (set-configuration)
          (nrs/start-server :port rpl-port :handler nrepl-handler))

        (options :help)
        (do
          (println summary)
          (System/exit 0))

        :else
        (do (println "Unknown options " options)
            (System/exit 1))))))
