(defproject panclus "0.1.0-SNAPSHOT"
  :description "Pangenome clustering exploration"
  :url "http://example.com/FIXME"
  :license {:name "MIT"
            :url "http://opensource.org/licenses/MIT"}
  :dependencies [[org.clojure/clojure "1.8.0"]
                 [lambda-ml           "0.1.0"]
                 [hiccup "1.0.5"]
                 [expresso "0.2.2"]
                 [aysylu/loom "0.5.4"]
                 [loom-gorilla "0.1.0"]
                 [gorilla-renderable "2.0.0"]
                 [gorilla-plot "0.1.4"]

                 [gyptis "0.2.2"]

                 [org.clojure/data.csv "0.1.2"]
                 [net.apribase/clj-dns "0.1.0"]
                 [aerial.fs "1.1.5"]
                 [aerial.utils "1.2.0"]
                 [aerial.bio.utils "2.0.0"]]

  :main ^:skip-aot gorilla-test.core
  :target-path "target/%s"

  :plugins [[lein-gorilla "0.4.0" :exclusions [org.clojure/clojure]]
            ^:replace
            [refactor-nrepl    "2.2.0"]
            [cider/cider-nrepl "0.12.0"]]

  :profiles {:uberjar {:aot :all}}

  :repositories [["lclrepo" "file:lclrepo"]]
  :source-paths ["src"]
  )

