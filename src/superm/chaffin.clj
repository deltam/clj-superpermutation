(ns superm.core
  "Find minimal superpermutation with Chaffin Method"
  (:require [clojure.math.combinatorics :as cmb]
            [clojure.set :as cs]
            [clojure.core.async :as ac]))

(def +all-perm-set+ (map #(set (cmb/permutations (range 1 (inc %))))
                         (range)))

(defn perm-set [n]
  (nth +all-perm-set+ n))

(defn count-perms [n]
  (cmb/count-permutations (range 1 (inc n))))


(defn edge-cost
  "Calculate num of adding digits that concat p1->p2 "
  [p1 p2]
  (loop [c 1]
    (let [p1-suffix (drop c p1)]
      (if (every? (fn [[d1 d2]] (= d1 d2))
                  (map vector p1-suffix p2))
        c
        (recur (inc c))))))

(defn rotate-origin
  "Return permutation of rotate `perm` that start with 1"
  [perm]
  (let [s (concat perm perm)]
    (take (count perm) (drop-while #(not= % 1) s))))

(defn rot= [p1 p2 level]
  (let [n (count p1)
        f (fn [p] (filter #(<= % (- n level)) p))]
    (= (rotate-origin (f p1)) (rotate-origin (f p2)))))

(defn connect? [s d c]
  (let [n (count s)]
    (loop [i 0]
      (if (rot= s d i)
        (<= c (inc i))
        (if (< i (- n 3))
          (recur (inc i))
          (< c n))))))

(defn raw-perm-graph [n]
  (let [pset (perm-set n)
        es (for [s pset, d pset
                 :let [c (edge-cost s d)]
                 :when (and (not= s d)
                            (< c n)
                            (connect? s d c))]
          [[s d] c])]
    {:nodes pset
     :edges (apply hash-map
                   (apply concat es))}))

(def +all-perm-graph+ (map raw-perm-graph (range)))

(defn perm-graph [n]
  (nth +all-perm-graph+ n))


(defn gen-perm-seq [n]
  {:n n
   :rest (perm-set n)
   :seq []
   :waste 0})

(defn init-perm-seq
  ([n init]
   (let [ps (gen-perm-seq n)]
     (-> ps
         (update :rest disj init)
         (update :seq conj {:perm init, :cost n, :waste '()}))))
  ([n] (init-perm-seq n (vec (range 1 (inc n))))))

(defn last-perm [perm-seq]
  (:perm (last (:seq perm-seq))))

(defn tail-seq [perm-seq perm c]
  (let [tail (concat (drop 1 (last-perm perm-seq)) (take-last c perm))]
    (map #(->> tail
               (drop %)
               (take (:n perm-seq)))
         (range c))))

(defn waste-conj [perm-seq perm c]
  (let [tails (tail-seq perm-seq perm c)]
    (map last
         (filter #(not (contains? (:rest perm-seq) %))
                 tails))))

(defn conj-perm [perm-seq perm c]
  (let [w (waste-conj perm-seq perm c)
        all (perm-set (:n perm-seq))
        tail-perms (filter #(contains? all %) (tail-seq perm-seq perm c))]
    (-> perm-seq
        (update :seq conj {:perm perm, :cost c, :waste w})
        (update :rest cs/difference (set tail-perms))
        (update :waste + (count w)))))

(defn rank [perm-seq]
  (let [all (count-perms (:n perm-seq))]
    (- all (count (:rest perm-seq)))))

(defn find-next-perm [perm-seq]
  (let [edges (:edges (perm-graph (:n perm-seq)))
        lp (last-perm perm-seq)
        nexts (map #(vector lp %) (:rest perm-seq))]
    (select-keys edges nexts)))



(defn short-branch? [cur-perm max-perm table]
  (let [w (- (:waste max-perm) (:waste cur-perm))
        diff (- (rank max-perm) (rank cur-perm))]
    (and (<= 0 w) (< w (count table))
         (<= (nth table w) diff))))

(defn chaffin-search
  "Search for the longest perm-seq with waste digits below waste-limit"
  [prefix-perm max-perm waste-limit table]
  (let [branches (->> (find-next-perm prefix-perm)
                      (map (fn [[[_ d] c]] (conj-perm prefix-perm d c)))
                      (filter (fn [{w :waste}] (<= w waste-limit))))]
    (if (empty? branches)
      (if (< (rank max-perm) (rank prefix-perm))
        prefix-perm
        max-perm)
      (reduce (fn [mp p] (if (short-branch? p mp table)
                           mp
                           (chaffin-search p mp waste-limit table)))
              max-perm
              branches))))

(defn chaffin-search-same-length
  [prefix-perm max-perm sames]
  (let [branches (->> (find-next-perm prefix-perm)
                      (map (fn [[[_ d] c]] (conj-perm prefix-perm d c)))
                      (filter (fn [{w :waste}] (<= w (:waste max-perm)))))]
    (if (empty? branches)
      (if (= (rank max-perm) (rank prefix-perm))
        (conj sames prefix-perm)
        sames)
      (reduce (fn [sm p] (chaffin-search-same-length p max-perm sm))
              sames
              branches))))

(defn chaffin-seq [n]
  (let [ps (init-perm-seq n)
        cs (->> (iterate (fn [[_ w table]]
                           (let [mp (chaffin-search ps ps w table)]
                             [mp (inc w) (conj table (rank mp))]))
                         [nil 0 []])
                (rest)
                (map first))
        all (count-perms n)
        lim (take-while #(< (rank %) all)
                        cs)]
    (lazy-cat lim (take 1 (drop (count lim) cs)))))

(def +all-chaffin-table+ (let [all (map #(map rank (chaffin-seq %)) (range))]
                           (concat [[0] [1] [2] [3 6] [4 8 12 14 18 20 24]]
                                   (drop 5 all))))

(defn chaffin
  "Chaffin Method https://github.com/superpermutators/superperm/wiki/Chaffin-method"
  ([n waste-limit]
   (let [table (take waste-limit (nth +all-chaffin-table+ n))]
     (chaffin n waste-limit table)))
  ([n waste-limit table]
   (let [ps (init-perm-seq n)]
     (chaffin-search ps ps waste-limit table))))



(defn perms->digits [{perms :seq}]
  (reduce (fn [dig {p :perm, c :cost}]
            (concat dig (take-last c p)))
          (:perm (first perms))
          (rest perms)))

(defn perms->str [perm-seq]
  (apply str (perms->digits perm-seq)))

(defn str->perms [s n]
  (let [is (map #(Integer/parseInt (str %)) s)
        all (perm-set n)
        perms (filter #(contains? all %)
                      (take-while #(= (count %) n)
                                  (map #(take n (drop % is)) (range))))
        ps (init-perm-seq n (first perms))
        g (perm-graph n)]
    (loop [ret ps, p (first perms), p2 (second perms), r (drop 2 perms)]
      (if (nil? p2)
        ret
        (let [c ((:edges  g) [p p2])]
          (recur (conj-perm ret p2 c) p2 (first r) (rest r)))))))


(defn print-chaffin-table [n]
  (let [cs (map #(hash-map :waste (:waste %)
                           :max (rank %)
                           :perm (perms->str %))
                (chaffin-seq n))]
    (printf "n=%d\n" n)
    (doseq [{w :waste, m :max, p :perm} cs]
      (printf "%d\t%d\t%s\n" w m p))))

(defn print-subs
  "Print all substring of length n and mark `x` if not permutation"
  [ps]
  (let [n (:n ps)
        ds (perms->digits ps)]
    (loop [p (take n ds), r (rest ds), indent "", s (perm-set n)]
      (when (= n (count p))
        (print (str indent (apply str p)))
        (if (not (contains? s p))
          (print " x"))
        (println)
        (recur (take n r) (rest r) (str indent " ") (disj s p))))))



(comment

(perms->str (chaffin 3 1))
; "123121321"

(time (def ans (chaffin 4 6)))
; "Elapsed time: 71728.029332 msecs"
(perms->digits ans)
; "123412314231243121342132413214321"
(rank ans)
; 24


(print-chaffin-table 4)
; n=4
; 0	4	1234123
; 1	8	123412314231
; 2	12	12341231423124312
; 3	14	12341231423124321432
; 4	18	1234123142312432143241324
; 5	20	1234132413421341231423124312
; 6	24	123412314231243121342132413214321

;; https://github.com/superpermutators/superperm/blob/master/ChaffinMethodResults/ChaffinMethodMaxPerms_4.txt

)
