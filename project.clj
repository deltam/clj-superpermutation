(defproject clj-superpermutation "0.1.0-SNAPSHOT"
  :description "Find minimal super permutation"
  :url "https://github.com/deltam/clj-superpermutation"
  :license {:name "MIT License"
            :url "none"
            :year 2019
            :key "mit"}
  :dependencies [[org.clojure/clojure "1.10.0"]
                 [org.clojure/math.combinatorics "0.1.5"]
                 [org.clojure/core.async "0.4.490"]]
  :repl-options {:init-ns superperm.core})
