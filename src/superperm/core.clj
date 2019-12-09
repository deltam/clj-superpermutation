(ns superperm.core
  "Find minimal super permutations"
  (:require [clojure.math.combinatorics :as cmb]
            [clojure.set :as cs]))

(def +all-perm-set+ (map #(set (cmb/permutations (range 1 (inc %))))
                         (range)))

(defn perm-set [n]
  (nth +all-perm-set+ n))

(defn count-perms [n]
  (cmb/count-permutations (range 1 (inc n))))





(def ^:const perm-n 4)
(def upper-limit
  "n! + (n-1)! + (n-2)! + (n-3)! + n - 3
  https://www.gregegan.net/SCIENCE/Superpermutations/Superpermutations.html#WILLIAMS"
  (let [fn-3 (apply * (range 1 (- perm-n 2)))
        fn-2 (* fn-3 (- perm-n 2))
        fn-1 (* fn-2 (dec perm-n))
        fn (* fn-1 perm-n)]
    (+ fn fn-1 fn-2 fn-3 perm-n -3)))

(def pset (perm-set perm-n))

(defn cost [p1 p2]
  (first
   (filter #(= (drop % p1) (drop-last % p2))
           (range 1 (inc perm-n)))))

(def edges-cost (let [es (for [x pset, y pset :when (not= x y)]
                           [(cost x y) x y])]
                  (zipmap (map rest es) (map first es))))

(defn cost-under? [vs]
  (every? #(< % perm-n)
          (map cost vs (rest vs))))

(defn vec->spm [vs]
  (reduce (fn [spm p]
            (let [c (cost (take-last perm-n spm) p)]
              (concat spm (take-last c p))))
          (first vs)
          (rest vs)))

(def all-spm (let [st (range 1 (inc perm-n))]
               (->> (cmb/permutations (disj pset st))
                    (filter cost-under?)
                    (map #(vec->spm (concat [st] %))))))

(defn min-spm []
  (let [cut (filter #(<= (count %) upper-limit) all-spm)
        sorted (sort #(< (count %1) (count %2)) cut)
        mn (first sorted)]
    (take-while #(= (count %) (count mn)) sorted)))




(defn cost? [p1 p2 c]
  (= (drop c p1) (drop-last c p2)))

(defn gen-edges [ps cost]
  (zipmap ps
          (map (fn [p] (set (filter #(cost? p % cost) ps)))
               ps)))

(def edges (let [cs (range 1 perm-n)]
             (zipmap cs
                     (map #(gen-edges pset %)
                          cs))))

(defn next-perms [p ps]
  (mapcat (fn [[c es]]
            (map #(vector c %) (filter ps (es p))))
          edges))

(defn spm->perms [spm]
  (->> (range)
       (map #(take perm-n (drop % spm)))
       (take-while #(= (count %) perm-n))
       (set)))

(defn find-superperm
  ([spm_ mn c p]
   (let [spm (concat spm_ (take-last c p))
         rest-ps (cs/difference pset (spm->perms spm))]
     (cond
       (< (count mn) (count spm)) mn
       (empty? rest-ps) spm
       :else (reduce (fn [m [c p]]
                       (let [m2 (find-superperm spm m c p)]
                         (if (< (count m2) (count m))
                           m2
                           m)))
                     mn
                     (next-perms (take-last perm-n spm) rest-ps)))))
  ([] (let [st (range 1 (inc perm-n))]
        (find-superperm [] (range 100) perm-n st))))


(declare find-chaffin)
(def chaffin-table (map #(:rank (find-chaffin %)) (range)))

(defn find-chaffin
  ([spm_ ps_ waste_ mx-waste mx c p]
   (let [spm (concat spm_ (take-last c p))
         add-ps (spm->perms (take-last (+ c perm-n -1) spm))
         ps (cs/difference ps_ add-ps)
         waste (+ waste_ (count (cs/difference add-ps ps_)))
         rank (- (count pset) (count ps))
         rw (- mx-waste waste)]
     (cond
       (< mx-waste waste)         {:spm spm_, :rank (- (count pset) (count ps_))}
       (< (+ rank (nth chaffin-table rw (count pset))) (:rank mx)) mx
       (empty? ps)                {:spm spm, :rank rank}
       :else (reduce (fn [m [c p]]
                       (let [m2 (find-chaffin spm ps waste mx-waste m c p)]
                         (if (< (:rank m) (:rank m2))
                           m2
                           m)))
                     mx
                     (next-perms (take-last perm-n spm) ps)))))
  ([mx-waste] (let [st (range 1 (inc perm-n))]
                (find-chaffin [] pset 0 mx-waste {:rank 0} perm-n st))))
