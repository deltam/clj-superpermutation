(ns superperm.core
  "Find minimal super permutations"
  (:require [superperm.tree :refer :all]
            [clojure.math.combinatorics :as cmb]
            [clojure.set :as cs]))

(def +all-perm-set+ (map #(set (cmb/permutations (range 1 (inc %))))
                         (range)))

(defn perm-set [n]
  (nth +all-perm-set+ n))

(defn count-perms [n]
  (cmb/count-permutations (range 1 (inc n))))


;; perm config

(defn raw-cost [p1 p2]
  (let [n (count p1)]
    (first
     (filter #(= (drop % p1) (drop-last % p2))
             (range 1 (inc n))))))

(defn upper-limit [n]
  "n! + (n-1)! + (n-2)! + (n-3)! + n - 3
  https://www.gregegan.net/SCIENCE/Superpermutations/Superpermutations.html#WILLIAMS"
  (let [fn-3 (apply * (range 1 (- n 2)))
        fn-2 (* fn-3 (- n 2))
        fn-1 (* fn-2 (dec n))
        fn (* fn-1 n)]
    (+ fn fn-1 fn-2 fn-3 n -3)))

(defn gen-config [n]
  {:n n
   :pset (perm-set n)
   :digits (range 1 (inc n))
   :upper-limit (upper-limit n)
   :costf (memoize raw-cost)})

(def ^:dynamic *config* (atom nil))

(defn reset-config! [n] (reset! *config* (gen-config n)))
(reset-config! 3)

(defn cfg [k]
  (if (= k :pset-count)
    (count (:pset @*config*))
    (@*config* k)))

(defn cost [a b] ((cfg :costf) a b))



;; superperm tree

(defn spm->perms [spm]
  (let [n (cfg :n)]
    (->> (range)
         (map #(take n (drop % spm)))
         (take-while #(= (count %) n)))))

(defn conj-perm [prefix p]
  (let [n (cfg :n)
        c (cost (take-last n prefix) p)
        cur (concat prefix (take-last c p))
        conj-ps (set (spm->perms (take-last (+ n -1 c) cur)))]
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

(defn spm-tree []
  (rep-tree gen-spm-branch {:spm (cfg :digits)
                            :rest (disj (cfg :pset) (cfg :digits))}))

(defn spm? [v] (empty? (:rest v)))

(defn find-min-spm [t]
  (prune-reduce-leaves (fn [mn v] (<= (count mn) (count (:spm v))))
                       (fn [mn v] (if (spm? v)
                                    (:spm v)
                                    mn))
                       (range (cfg :upper-limit))
                       t))

(defn conj-min [ms m]
  (filter #(<= (count %) (count m))
          (conj ms m)))

(defn find-min-spms [t]
  (first
   (prune-reduce-leaves (fn [[mns mc] v] (< mc (count (:spm v))))
                        (fn [[mns mc] v] (if (spm? v)
                                           [(conj-min mns (:spm v)) (count (:spm v))]
                                           [mns mc]))
                        [[] (cfg :upper-limit)]
                        t)))


;; utils

(defn tree->dot
  ([f g t]
   (let [cnt (atom 0)
         [edges nodes _]
         (reduce-tree (fn [[edges nodes ids bs] v]
                        (let [lbl (f v)
                              id (format "v%s_%d" (apply str (:spm v)) @cnt)
                              nd (if (empty? (:rest v))
                                   (format "%s [label=\"%s\"; shape=box; color=blue;];" id lbl)
                                   (format "%s [label=\"%s\";];" id lbl))]
                          (swap! cnt inc)
                          [(concat edges (map #(format "%s -> %s[label=\"%s\"];" id %1 (g v %2)) ids bs))
                           (conj nodes nd)
                           id
                           v]))
                      (fn [[edges nodes ids bs] [es ns id b]]
                        [(concat edges es)
                         (concat nodes ns)
                         (conj ids id)
                         (conj bs b)])
                      [[] [] [] []]
                      t)]
     (format "digraph tree{\n graph[rankdir=LR;];\n %s \n %s}"
             (clojure.string/join "\n" nodes)
             (clojure.string/join "\n" edges))))
  ([t] (tree->dot #(apply str (:spm %))
                  #(apply str (take-last (cfg :n) (:spm %2)))
                  t)))
