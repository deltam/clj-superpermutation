(ns superm.core
  (:require [clojure.math.combinatorics :as cmb]
            [clojure.core.async :as ac]))

(def +all-perm-set+ (map #(apply hash-set
                                 (cmb/permutations (range 1 (inc %))))
                         (range)))

(defn perms [n]
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

(defn cost-graph [n perms]
  (let [pairs (for [p1 perms, p2 perms :when (not= p1 p2)]
                [p1 p2])]
    {:nodes (apply hash-set perms)
     :edges (apply merge
                   (for [[src dest] pairs
                         :let [c (edge-cost src dest)]
                         :when (< c n)]
                     {[src dest] c}))}))

(defn perm-graph [n]
  (cost-graph n (perms n)))

(defn cost [g s d]
  ((:edges g) [s d]))

(defn print-graph [g]
  (doseq [[[p1 p2] c] (:edges g)]
    (printf "%s -- %d --> %s\n" p1 c p2)))


(defn gen-perm-seq [n]
  {:n n
   :reached #{}
   :seq []
   :waste 0})

(defn init-perm-seq
  ([n init]
   (let [ps (gen-perm-seq n)]
     (-> ps
         (update :reached conj init)
         (update :seq conj {:perm init, :cost 0}))))
  ([n] (init-perm-seq n (vec (range 1 (inc n))))))

(defn last-perm [perm-seq]
  (:perm (last (:seq perm-seq))))

(defn waste-conj [perm-seq perm cost]
  (let [n (:n perm-seq)
        all (perms n)
        lp (last-perm perm-seq)
        sub-perm (concat lp (take-last cost perm))]
    (count
     (filter #(or (not (contains? all %))
                  (contains? (:reached perm-seq) %))
             (map #(drop-last % (take-last n sub-perm))
                  (range cost))))))

(defn conj-perm [perm-seq perm c]
  (-> perm-seq
      (update :seq conj {:perm perm, :cost c})
      (update :reached conj perm)
      (update :waste + (waste-conj perm-seq perm c))))

(defn rank [perm-seq]
  (count (:reached perm-seq)))

(defn find-next-perm [perm-seq g]
  (let [dest-nodes (clojure.set/difference (perms (:n perm-seq)) (:reached perm-seq))
        lp (last-perm perm-seq)
        edges (map #(vector lp %) dest-nodes)]
    (select-keys (:edges g) edges)))

(defn max-rank [ps]
  (reduce (fn [mp p] (if (< (rank mp) (rank p))
                       p
                       mp))
          ps))


(defn chaffin-branch-seq
  "Return coll of perm-seq has waste under waste-limit"
  [prefix-perm g waste-limit]
  (let [branches (map (fn [[[_ d] c]] (conj-perm prefix-perm d c))
                      (find-next-perm prefix-perm g))
        cut (filter (fn [{w :waste}] (<= w waste-limit))
                    branches)]
    (if (empty? cut)
      [prefix-perm]
      (mapcat #(chaffin-branch-seq % g waste-limit) cut))))

(defn chaffin
  "Chaffin Method https://github.com/superpermutators/superperm/wiki/Chaffin-method"
  [n waste-limit]
  (let [g (perm-graph n)
        ps (init-perm-seq n)]
    (max-rank (chaffin-branch-seq ps g waste-limit))))


(defn perms->digits [{perms :seq}]
  (reduce (fn [dig {p :perm, c :cost}]
            (concat dig (take-last c p)))
          (:perm (first perms))
          (rest perms)))

(defn perms->str [perm-seq]
  (apply str (perms->digits perm-seq)))

(defn print-chaffin-table [n]
  (let [chaffin-seq (map (fn [w] (let [c (chaffin n w)]
                                   {:waste w,
                                    :max (rank c)
                                    :perm (perms->str c)}))
                         (range))
        all (count-perms n)]
    (printf "n=%d\n" n)
    (loop [{w :waste, m :max, p :perm} (first chaffin-seq), cs (rest chaffin-seq)]
      (printf "%d\t%d\t%s\n" w m p)
      (if (< m all)
        (recur (first cs) (rest cs))))))

(defn str->perms [s n]
  (let [is (map #(Integer/parseInt (str %)) s)
        all (perms n)
        perms (filter #(contains? all %)
                      (take-while #(= (count %) n)
                                  (map #(take n (drop % is)) (range))))
        ps (init-perm-seq n (first perms))
        g (perm-graph n)]
    (loop [ret ps, p (first perms), p2 (second perms), r (drop 2 perms)]
      (if (nil? p2)
        ret
        (let [c (cost g p p2)]
          (recur (conj-perm ret p2 c) p2 (first r) (rest r)))))))



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
