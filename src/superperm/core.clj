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
(def pset-count (count pset))
(def perm-digits (range 1 (inc perm-n)))

(defn raw-cost [p1 p2]
  (first
   (filter #(= (drop % p1) (drop-last % p2))
           (range 1 (inc perm-n)))))
(def cost (memoize raw-cost))



;; functional tree

(defn node [v bs] (vector v bs))

(defn rep-tree [f v]
  (node v (fn [] (map #(rep-tree f %)
                      (f v)))))

(defn tval [t] (first t))
(defn branch [t] ((second t)))

(defn reduce-tree [f g init t]
  (f (reduce #(g %1 (reduce-tree f g init %2))
             init
             (branch t))
     (tval t)))

(defn map-tree [f t]
  (reduce-tree (fn [bs v] (node (f v) (fn [] bs)))
               conj
               []
               t))

(defn take-while-tree [pred t]
  (if (pred (tval t))
    (node (tval t) (fn [] (filter not-empty
                                  (map #(take-while-tree pred %)
                                       (branch t)))))))

(defn reduce-leaves [f init t]
  (let [bs (branch t)]
    (if (empty? bs)
      (f init (tval t))
      (reduce #(reduce-leaves f %1 %2) init bs))))



;; superperm tree

(defn spm->perms [spm]
  (->> (range)
       (map #(take perm-n (drop % spm)))
       (take-while #(= (count %) perm-n))
       (set)))

(defn conj-perm [prefix p]
  (let [c (cost (take-last perm-n prefix) p)
        cur (concat prefix (take-last c p))
        conj-ps (spm->perms (take-last (+ perm-n -1 c) cur))]
    [cur conj-ps]))

(defn pre= [a b]
  (let [l (min (count a) (count b))]
    (= (take l a) (take l b))))

(defn distinct-prefix [spms]
  (->> spms
       (sort #(< (count (:spm %1)) (count (:spm %2))))
       (reduce (fn [acc s]
                 (if (some #(pre= (:spm s) (:spm %)) acc)
                   acc
                   (conj acc s)))
               [])))

(defn gen-spm-branch [{prefix :spm, rest-ps :rest}]
  (->> rest-ps
       (map #(let [[cur conj-ps] (conj-perm prefix %)]
               {:spm cur, :rest (cs/difference rest-ps conj-ps)}))
       (distinct-prefix)))

(def spm-tree
  (rep-tree gen-spm-branch {:spm perm-digits
                            :rest (disj pset perm-digits)}))


(defn prune-over-branch [t]
  (take-while-tree (fn [{spm :spm}] (<= (count spm) upper-limit))
                   t))

(defn min-count [a b]
  (min-key count b a))

(defn find-min-spm [t]
  (reduce-leaves (fn [mn {spm :spm, ps :rest}]
                   (if (empty? ps)
                     (min-count spm mn)
                     mn))
                 (range upper-limit)
                 t))

; (->> spm-tree (prune-over-branch) (find-min-spm))

(defn filter-min-count [vs]
  (let [mn (reduce min-count vs)]
    (filter #(<= (count %) (count mn)) vs)))

(defn find-min-spms [t]
  (distinct
   (reduce-leaves (fn [mns {spm :spm, ps :rest}]
                    (if (empty? ps)
                      (filter-min-count (conj mns spm))
                      mns))
                  []
                  t)))



;; chaffin method

(defn gen-chaffin-branch [{prefix :spm, waste :waste, rest-ps :rest}]
  (->> rest-ps
       (map #(let [[cur conj-ps] (conj-perm prefix %)
                   w (count (cs/difference conj-ps rest-ps))]
               {:spm cur
                :waste (+ waste w)
                :rest (cs/difference rest-ps conj-ps)}))
       (distinct-prefix)))

(def chaffin-tree
  (rep-tree gen-chaffin-branch {:spm perm-digits
                                :waste 0
                                :rest (disj pset perm-digits)}))

(defn take-while-waste [waste t]
  (take-while-tree (fn [{w :waste}] (<= w waste))
                   t))

(defn max-contain-perm [mv pv]
  (max-key first pv mv))

(defn find-max-contain-perm [t]
  (reduce-leaves (fn [mx {spm :spm, rest-ps :rest}]
                   (let [c (- pset-count (count rest-ps))]
                     (max-contain-perm mx [c spm])))
                 [0 []]
                 t))

(defn find-chaffin [w]
  (->> chaffin-tree
       (take-while-waste w)
       (find-max-contain-perm)))


;; utils

(defn print-tree [f t]
  (reduce-tree (fn [ac v] (conj [(f v)] ac))
               conj
               []
               t))

(defn prune [n t]
  (if (< 0 n)
    (node (tval t) (fn [] (filter not-empty
                                  (map #(prune (dec n) %) (branch t)))))))

(defn get-tree [t [k & ks]]
  (let [bs (branch t)]
    (cond
      (nil? k) t
      (empty? bs) []
      :else (get-tree (nth bs k) ks))))

(defn count-leaves [t]
  (reduce-leaves (fn [cnt _] (inc cnt))
                 0
                 t))
