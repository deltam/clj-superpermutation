(ns superperm.chaffin-test
  (:require [clojure.test :refer :all]
            [superperm.chaffin :refer :all]))

(deftest tail-seq-test
  (let [perm-seq (init-perm-seq 3)]
    (is (= '((2 3 1))
           (tail-seq perm-seq [2 3 1] 1)))
    (is (= '((2 3 1) (3 1 2))
           (tail-seq perm-seq [3 1 2] 2)))))

(deftest waste-conj-test
  (let [ps (-> (init-perm-seq 3)
               (conj-perm [2 3 1] 1)
               (conj-perm [3 1 2] 1))]
    (is (= '(1) (waste-conj ps [2 1 3] 2)))))

(deftest conj-perm-test
  (let [ps (init-perm-seq 3)]
    (is (= {:n 3
            :reached #{[1 2 3] [2 3 1]}
            :seq [{:perm [1 2 3], :cost 0, :waste '()}
                  {:perm [2 3 1], :cost 1, :waste '()}]
            :waste 0}
           (conj-perm ps [2 3 1] 1)))
    (is (= {:n 3
            :reached #{[1 2 3] [3 2 1]}
            :seq [{:perm [1 2 3], :cost 0, :waste '()}
                  {:perm [3 2 1], :cost 2, :waste '(2)}]
            :waste 1}
           (conj-perm ps [3 2 1] 2)))))
