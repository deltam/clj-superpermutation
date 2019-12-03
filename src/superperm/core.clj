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

(def ^:const perm-n 4)
(def pset (perm-set perm-n))
(def edges (let [cs (range 1 perm-n)]
             (zipmap cs
                     (map #(gen-edges pset %)
                          cs))))

(defn next-perms [p ps]
  (mapcat (fn [[c es]]
            (map #(vector c %) (filter ps (es p))))
          edges))

(defn waste-count [p1 p2 c ps]
  (let [p (concat p1 (take-last c p2))]
    (->> (map #(take-last perm-n (drop-last % p)) (range c))
         (filter (complement ps))
         (count))))

(defn find-superperm
  ([spm cost waste mx-waste mx ps]
   (cond
     (empty? ps) {:spm spm, :cost cost, :waste waste}
     :else (reduce (fn [m [c p]]
                     (let [w (+ waste (waste-count (last spm) p c ps))
                           m2 (if (< mx-waste w)
                                {:spm spm, :cost cost, :waste waste}
                                (find-superperm (conj spm p) (+ cost c) w mx-waste m (disj ps p)))]
                       (if (< (:cost m) (:cost m2))
                         m2
                         m)))
                   mx
                   (next-perms (last spm) ps))))
  ([mx-waste] (let [st (range 1 (inc perm-n))]
        (find-superperm [st] 0 0 mx-waste {:cost 0} (disj pset st)))))

(defn print-superperm [spm]
  (apply str
         (first
          (reduce (fn [[r lp] p]
                    (let [c (first (filter #(cost? lp p %) (range perm-n)))]
                      [(concat r (take-last c p)) p]))
                  [(first spm) (first spm)]
                  (rest spm)))))

;(find-superperm [[1 2 3]] 0 pset)
