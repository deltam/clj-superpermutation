(ns superperm.tree)

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

(defn prune-reduce-leaves [p f init t]
  (let [acc (atom init)]
    (->> t
         (take-while-tree #(not (p @acc %)))
         (reduce-leaves #(do (reset! acc (f %1 %2))
                             @acc)
                        init))))



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

(defn count-leaves [t]
  (reduce-leaves (fn [cnt _] (inc cnt))
                 0
                 t))
