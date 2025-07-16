(define-module (mittleff utils)
  #:use-module (srfi srfi-1)
  #:use-module (system foreign)
  #:use-module (system foreign-library)
  #:export (accuracy-from-prec))

;; Load the shared library
(define libutils
  (catch #t
    (lambda () (load-foreign-library "libutils" #:global? #t))
    (lambda (key . args)
      (error "Failed to load libutils:" args))))

;; Keep the compute-r1 function as it was (not moved to C library)
(define* (accuracy-from-prec prec)
  (let ((func (foreign-library-function libutils "accuracy_from_prec"
                  #:return-type double
                  #:arg-types (list int))))
    (func prec)))
