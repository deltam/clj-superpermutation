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





(def ^:const perm-n 3)
(def ^:const upper-limit
  "n! + (n-1)! + (n-2)! + (n-3)! + n - 3
  https://www.gregegan.net/SCIENCE/Superpermutations/Superpermutations.html#WILLIAMS"
  (let [fn-3 (apply * (range 1 (- perm-n 2)))
        fn-2 (* fn-3 (- perm-n 2))
        fn-1 (* fn-2 (dec perm-n))
        fn (* fn-1 perm-n)]
    (+ fn fn-1 fn-2 fn-3 perm-n -3)))

(def pset (perm-set perm-n))

(defn raw-cost [p1 p2]
  (first
   (filter #(= (drop % p1) (drop-last % p2))
           (range 1 (inc perm-n)))))
(def cost (memoize raw-cost))



;; functional tree

(defn rep-tree [f v]
  [v (fn [] (map #(rep-tree f %)
                 (f v)))])

(defn tval [t] (first t))
(defn branch [t] ((second t)))

(declare reduce-tree-bs)
(defn reduce-tree [f g init t]
  (f (reduce-tree-bs f g init (branch t))
     (tval t)))
(defn reduce-tree-bs [f g init bs]
  (if (empty? bs)
    init
    (g (reduce-tree f g init (first bs))
       (reduce-tree-bs f g init (rest bs)))))

(defn take-while-tree [pred t]
  (cond
    (empty? t) t
    (pred (tval t)) []
    :else [(tval t) (fn [] (filter not-empty
                                   (map #(take-while-tree pred %) (branch t))))]))



;; superperm tree

(defn gen-spm-branch [[prefix rest-ps]]
  (map #(let [c (cost (take-last perm-n prefix) %)]
          [(concat prefix (take-last c %)) (disj rest-ps %)])
       rest-ps))

(def spm-tree (let [st (range 1 (inc perm-n))]
                (rep-tree gen-spm-branch [st (disj pset st)])))

(defn min-count [a b]
  (if (< (count a) (count b))
    a
    b))

(defn find-min-spm [t]
  (->> t
       (take-while-tree (fn [[spm _]] (< upper-limit (count spm))))
       (reduce-tree (fn [mn [spm ps]] (if (empty? ps)
                                        (min-count spm mn)
                                        mn))
                    min-count
                    (range upper-limit))))


(defn filter-min-count [vs]
  (let [mn (reduce min-count (first vs) (rest vs))]
    (filter #(<= (count %) (count mn)) vs)))

(defn find-min-spms [t]
  (->> t
       (take-while-tree (fn [[spm _]] (< upper-limit (count spm))))
       (reduce-tree (fn [ms [spm ps]] (if (empty? ps)
                                        (conj ms spm)
                                        ms))
                    (fn [ms spms] (filter-min-count (concat ms spms)))
                    [])))



;; utils

(defn prune [n t]
  (if (<= n 0)
    []
    [(tval t) (fn [] (filter not-empty
                             (map #(prune (dec n) %) (branch t))))]))

(defn get-tree [t [k & ks]]
  (let [bs (branch t)]
    (cond
      (nil? k) t
      (empty? bs) []
      :else (get-tree (nth bs k) ks))))

(defn count-leaf [t]
  (reduce-tree (fn [cnt [spm ps]] (if (empty? ps)
                                    (inc cnt)
                                    cnt))
               +
               0
               t))



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
