(ns superperm.core
  "Find minimal super permutations"
  (:require [clojure.math.combinatorics :as cmb]))

(def +all-perm-set+ (map #(set (cmb/permutations (range 1 (inc %))))
                         (range)))

(defn perm-set [n]
  (nth +all-perm-set+ n))

(defn count-perms [n]
  (cmb/count-permutations (range 1 (inc n))))
