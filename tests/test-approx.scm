(define-module (tests test-approx)
  #:use-module (srfi srfi-64)
  #:export (test-approx))

;; From https://numpy.org/doc/stable/reference/generated/numpy.isclose.html
(define* (numpy-isclose a b #:key (atol 1e-8) (rtol 1e-5))
  (<= (magnitude (- a b)) (+ atol (* rtol (max (magnitude a) (magnitude b))))))

(define* (test-approx msg expected computed #:key (atol 1e-8) (rtol 1e-5))
  (test-assert
      (format #f "~a\n-> Expected: ~a\n-> Computed: ~a\n" msg expected computed)
    (numpy-isclose computed expected)))


;; (define (relerr1 expected computed)
;;   (/ (magnitude (- expected computed)) (magnitude expected)))

;; ;; https://stackoverflow.com/a/34233736/3543586
;; (define* (test-approx1 msg expected computed #:optional (eps 1e-15))
;;   (let ((ulp-n 4)
;;         (diff (magnitude (- expected computed)))
;;         (mag  (* 0.5 (+ (magnitude expected) (magnitude computed)))))
;;     (test-assert (format #f "~a\n-> Expected: ~a\n-> Computed: ~a\n-> RelErr: ~a\n"
;;                          msg
;;                          expected
;;                          computed
;;                          (relerr expected computed))
;;       (or (< (relerr expected computed) 1e-12) (<= diff (* mag eps (ash 1 ulp-n)))))))



