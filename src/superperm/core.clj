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

(defn reduce-leaves [f init t]
  (let [bs (branch t)]
    (if (empty? bs)
      (f init (tval t))
      (reduce #(reduce-leaves f %1 %2) init bs))))

(defn map-tree [f t]
  (node (f (tval t))
        (fn [] (map #(map-tree f %) (branch t)))))

(defn take-while-tree [pred t]
  (if (pred (tval t))
    (node (tval t) (fn [] (filter not-empty
                                  (map #(take-while-tree pred %)
                                       (branch t)))))))




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

(defn gen-spm-tree []
  (rep-tree gen-spm-branch {:spm (cfg :digits)
                            :rest (disj (cfg :pset) (cfg :digits))}))

(defn prune-over-branch [t]
  (take-while-tree (fn [{spm :spm}] (<= (count spm) (cfg :upper-limit)))
                   t))

(defn min-count [a b]
  (min-key count b a))

(defn find-shortest-spm [t]
  (reduce-leaves (fn [mn {spm :spm, ps :rest}]
                   (if (empty? ps)
                     (min-count spm mn)
                     mn))
                 (range (cfg :upper-limit))
                 t))

(defn filter-min-count [vs]
  (let [mn (reduce min-count vs)]
    (filter #(<= (count %) (count mn)) vs)))

(defn find-shortest-spms [t]
  (distinct
   (reduce-leaves (fn [mns {spm :spm, ps :rest}]
                    (if (empty? ps)
                      (filter-min-count (conj mns spm))
                      mns))
                  []
                  t)))



;; chaffin method

(defn count-waste [spm]
  (let [ps (spm->perms spm)
        r (reduce disj (cfg :pset) ps)]
    (- (count ps)
       (- (cfg :pset-count) (count r)))))

(defn gen-chaffin-tree []
  (map-tree #(assoc % :waste (count-waste (:spm %)))
            (gen-spm-tree)))

(defn take-while-waste [waste t]
  (take-while-tree (fn [{w :waste}] (<= w waste))
                   t))

(defn max-contain-perm [mv pv]
  (max-key first pv mv))

(defn find-max-contain-perm [t]
  (reduce-leaves (fn [mx {spm :spm, rest-ps :rest}]
                   (let [c (- (cfg :pset-count) (count rest-ps))]
                     (max-contain-perm mx [c spm])))
                 [0 []]
                 t))

(defn find-chaffin [w]
  (->> (gen-chaffin-tree)
       (prune-over-branch)
       (take-while-waste w)
       (find-max-contain-perm)))

(defn print-chaffin-table []
  (let [cs (->> (range)
                (map #(concat [%] (find-chaffin %))))]
    (loop [[i m p] (first cs), r (rest cs)]
      (printf "%d\t%d\t%s\n" i m (apply str p))
      (if (< m (cfg :pset-count))
        (recur (first r) (rest r))))))




;; utils

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

(defn print-tree [f t]
  (reduce-tree (fn [ac v] (conj [(f v)] ac))
               conj
               []
               t))

(defn tree->dot [f t]
  (let [cnt (atom 0)
        [edges nodes _]
        (reduce-tree (fn [[edges nodes bs] v]
                       (let [lbl (f v)
                             id (format "v%s_%d" (apply str (:spm v)) @cnt)
                             nd (if (empty? (:rest v))
                                  (format "%s [label=\"%s\"; color=red; shape=box];" id lbl)
                                  (format "%s [label=\"%s\";];" id lbl))]
                         (swap! cnt inc)
                         [(concat edges (map #(format "%s -> %s;" id %) bs))
                          (conj nodes nd)
                          id]))
                     (fn [[edges nodes bs] [es ns s d]]
                       [(concat edges es)
                        (concat nodes ns)
                        (conj bs s)])
                     [[] [] []]
                     t)]
    (format "digraph tree{\n graph[rankdir=LR;];\n %s \n %s}"
            (clojure.string/join "\n" nodes)
            (clojure.string/join "\n" edges))))

;(->> spm-tree (prune-over-branch) (tree->dot #(apply str (:spm %))) (spit "tree3.dot"))

(defn count-leaves [t]
  (reduce-leaves (fn [cnt _] (inc cnt))
                 0
                 t))
