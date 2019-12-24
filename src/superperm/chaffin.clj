(ns superperm.chaffin
  "Implementation of Chaffin Method
  https://github.com/superpermutators/superperm/wiki/Chaffin-method"
  (:use [superperm.core :refer (cfg map-tree spm-tree spm->perms take-while-tree prune-reduce-leaves tval)]))


(defn count-waste [spm]
  (let [ps (spm->perms spm)
        r (reduce disj (cfg :pset) ps)]
    (- (count ps)
       (- (cfg :pset-count) (count r)))))

(defn chaffin-tree []
  (map-tree #(assoc %
                    :waste (count-waste (:spm %))
                    :rank (- (cfg :pset-count) (count (:rest %))))
            (spm-tree)))

(defn take-while-waste [waste t]
  (take-while-tree (fn [{w :waste}] (<= w waste))
                   t))

(defn find-chaffin [w tbl t]
  (->> t
       (take-while-waste w)
       (prune-reduce-leaves (fn [mx v] (let [dw (- (:waste mx) (:waste v))
                                             dr (- (:rank mx) (:rank v))]
                                         (<= (nth tbl dw (cfg :pset-count)) dr)))
                            (fn [mx v] (max-key :rank v mx))
                            (tval t))
       (#(select-keys % [:rank :spm :waste]))))

(defn chaffin-seq [t]
  (let [s (->> (iterate (fn [[_ tbl w]]
                          (let [c (find-chaffin w tbl t)]
                            [c (conj tbl (:rank c)) (inc w)]))
                        [nil [] 0])
               (rest)
               (map first))
        [tw dw] (split-with #(< (:rank %) (cfg :pset-count)) s)]
    (lazy-cat tw [(first dw)])))

(defn chaffin-table [t] (map :rank (chaffin-seq t)))

(defn print-chaffin-table []
  (doseq [{w :waste, r :rank, p :spm} (chaffin-seq (chaffin-tree))]
    (printf "%d\t%d\t%s\n" w r (apply str p))))






(comment

(time (superperm.chaffin/print-chaffin-table))
;> 0	4	1234123
;> 1	8	123412314231
;> 2	12	12341231423124312
;> 3	14	12341231423124321432
;> 4	18	1234123142312432143241324
;> 5	20	1234123142312431213421324132
;> 6	24	123412314231243121342132413214321
;> "Elapsed time: 1496.26679 msecs"

;; https://github.com/superpermutators/superperm/blob/master/ChaffinMethodResults/ChaffinMethodMaxPerms_4.txt

)
