(defproject panclus "0.1.0-SNAPSHOT"
  :description "Pangenome clustering exploration"
  :url "http://example.com/FIXME"
  :license {:name "MIT"
            :url "http://opensource.org/licenses/MIT"}
  :dependencies [[org.clojure/clojure "1.8.0"]
                 [org.clojure/tools.reader "1.0.0-beta3"]
                 [org.clojure/tools.nrepl  "0.2.13"] ; Explicit nREPL
                 [org.clojure/tools.cli    "0.3.3"] ; cmd line arg processing

                 [lambda-ml "0.1.0"]
                 [hiccup "1.0.5"]
                 [expresso "0.2.2"]
                 [aysylu/loom "0.5.4"]
                 [kixi/stats "0.4.0"]

                 [vgl "0.1.0"]

                 [org.clojure/data.csv "0.1.2"]
                 [net.apribase/clj-dns "0.1.0"]
                 [aerial.fs "1.1.5"]
                 [aerial.utils "1.2.0"]
                 [aerial.bio.utils "2.0.0"]]

  :target-path "target/%s"

  :plugins [^:replace
            [cider/cider-nrepl "0.17.0"]
            [refactor-nrepl    "2.4.0-SNAPSHOT"]]

  :profiles {:uberjar {:aot :all}}
  :main panclus.core

  :repositories [["lclrepo" "file:lclrepo"]]
  :source-paths ["src"]
  )

