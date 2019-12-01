(ns superperm.core
  "Find minimal super permutations"
  (:require [clojure.math.combinatorics :as cmb]))

(def +all-perm-set+ (map #(set (cmb/permutations (range 1 (inc %))))
                         (range)))

(defn perm-set [n]
  (nth +all-perm-set+ n))

(defn count-perms [n]
  (cmb/count-permutations (range 1 (inc n))))





(defn cost? [p1 p2 c]
  (= (drop c p1) (drop-last c p2)))

(defn gen-edges [ps cost]
  (zipmap ps
          (map (fn [p] (set (filter #(cost? p % cost) ps)))
               ps)))

(def ^:const perm-n 3)
(def pset (perm-set perm-n))
(def edges (let [cs (range 1 perm-n)]
             (zipmap cs
                     (map #(gen-edges pset %)
                          cs))))

(defn next-perms [p ps]
  (mapcat (fn [[c es]]
            (map #(vector c %) (filter ps (es p))))
          edges))

(defn find-superperm
  ([spm cost mn ps]
   (cond
     (< (:cost mn) cost) mn
     (empty? ps)         {:spm spm, :cost cost}
     :else (reduce (fn [m [c p]]
                     (let [m2 (find-superperm (conj spm p) (+ cost c) m (disj ps p))]
                       (if (< (:cost m2) (:cost m))
                         m2
                         m)))
                   mn
                   (next-perms (last spm) ps))))
;  (min-spm
;   (map (fn [[c p]] (find-superperm (conj spm p) (+ cost c) mn rest-ps))
;        es))
  ([st] (find-superperm [st] 0 {:cost Integer/MAX_VALUE} (disj pset st))))

(defn print-superperm [spm]
  (apply str
         (first
          (reduce (fn [[r lp] p]
                    (cond
                      (cost? lp p 1) [(concat r (take-last 1 p)) p]
                      (cost? lp p 2) [(concat r (take-last 2 p)) p]))
                  [(first spm) (first spm)]
                  (rest spm)))))

;(find-superperm [[1 2 3]] 0 pset)
