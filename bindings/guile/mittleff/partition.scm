(define-module (mittleff partition)
  ;;#:use-module (mittleff constants)
  #:use-module (srfi srfi-1)
  #:use-module (system foreign)
  #:use-module (system foreign-library)
  #:export (region?))

;; Load the shared library
;; (define libpartition
;;   (load-foreign-library "libpartition" #:global? #t))
(define libpartition
  (catch #t
    (lambda () (load-foreign-library "libpartition" #:global? #t))
    (lambda (key . args)
      (error "Failed to load libpartition:" args))))

(define* (region? z a taylor-radius prec)
  (let ((fn (foreign-library-function libpartition "region"
                            #:return-type unsigned-short
                            #:arg-types (list complex-double double double long))))
    (fn z a taylor-radius prec)))

;; Keep the compute-r1 function as it was (not moved to C library)
;; (define* (compute-r1 a #:key (prec 53))
;;   ;; (let* ((acc (expt 2 (- prec)))
;;   ;;        (pi 3.141592653589793238462643383279502884197169399375105820974944592307816406286)
;;   ;;       (c0 (/ (expt 1.3 (- 1 a)) (* pi (sin (* pi a))))))
;;   ;;   (expt (* -2 (log (/ acc c0))) a))
;;   (let ((func (foreign-library-function libpartition "compute_r1"
;;                   #:return-type double
;;                   #:arg-types (list double int))))
;;     (func a prec)))

;; ;; Factory function for region predicates
;; (define (make-region-predicate region-name)
;;   (lambda* (z a r0 #:key (prec 53))
;;     (let ((func (foreign-library-function libpartition region-name
;;                   #:return-type unsigned-short
;;                   #:arg-types (list complex-double double double long))))
;;       (= 1 (func z a r0 prec)))))

;; ;; Define all region predicates
;; (define in-region-G0? (make-region-predicate "in_region_G0"))
;; (define in-region-G1? (make-region-predicate "in_region_G1"))
;; (define in-region-G2? (make-region-predicate "in_region_G2"))
;; (define in-region-G3? (make-region-predicate "in_region_G3"))
;; (define in-region-G4? (make-region-predicate "in_region_G4"))
;; (define in-region-G5? (make-region-predicate "in_region_G5"))
;; (define in-region-G6? (make-region-predicate "in_region_G6"))

;; (define* (region z a r0 #:key (prec 53))
;;   "Determine which region the complex number z falls into for given parameter a.

;;    @var{z} The complex number to evaluate
;;    @var{a} The fractional order parameter (0 < a < 1)
;;    @var{prec} The precision in bits (default: 53)
;;    @var{radius} The boundary radius (default: *default-radius*)
;;    @return An integer from 0-6 representing the region"
;;   (cond
;;      ((in-region-G0? z a r0 #:prec prec) 0)
;;      (else
;;       (let ((r1 (compute-r1 a #:prec prec)))
;;         (cond
;;          ((>= (magnitude z) r1)
;;           ;; Check regions 1-4 in order of likely frequency
;;           (cond ((in-region-G1? z a r0 #:prec prec) 1)
;;                 ((in-region-G2? z a r0 #:prec prec) 2)
;;                 ((in-region-G3? z a r0 #:prec prec) 3)
;;                 ((in-region-G4? z a r0 #:prec prec) 4)
;;                 (else (error "Failed to classify z in regions 1-4"))))         
;;          (else
;;           ;; Check regions 5-6
;;           (cond ((in-region-G5? z a r0 #:prec prec) 5)
;;                 ((in-region-G6? z a r0 #:prec prec) 6)
;;                 (else (error "Failed to classify z in regions 5-6")))))))))

;; ;; Wrapper functions that call the C shared library
;; (define* (in-region-G0? z a #:key (prec 53))
;;   "Check if complex number z is in region G0"
;;   (let ((in-region-G0-c (foreign-library-function libpartition "in_region_G0"
;;                             #:return-type unsigned-short
;;                             #:arg-types (list complex-double double long))))
;;     (= 1 (in-region-G0-c z a prec))))

;; (define* (in-region-G1? z a #:key (prec 53))
;;   "Check if complex number z is in region G1"
;;   (let ((in-region-G1-c (foreign-library-function libpartition "in_region_G1"
;;                             #:return-type unsigned-short
;;                             #:arg-types (list complex-double double long))))
;;     (= 1 (in-region-G1-c z a prec))))

;; (define* (in-region-G2? z a #:key (prec 53))
;;   "Check if complex number z is in region G2"
;;   (let ((in-region-G2-c (foreign-library-function libpartition "in_region_G2"
;;                             #:return-type unsigned-short
;;                             #:arg-types (list complex-double double long))))
;;     (= 1 (in-region-G2-c z a prec))))

;; (define* (in-region-G3? z a #:key (prec 53))
;;   "Check if complex number z is in region G3"
;;   (let ((in-region-G3-c (foreign-library-function libpartition "in_region_G3"
;;                             #:return-type unsigned-short
;;                             #:arg-types (list complex-double double long))))
;;     (= 1 (in-region-G3-c z a prec))))

;; (define* (in-region-G4? z a #:key (prec 53))
;;   "Check if complex number z is in region G4"
;;   (let ((in-region-G4-c (foreign-library-function libpartition "in_region_G4"
;;                             #:return-type unsigned-short
;;                             #:arg-types (list complex-double double long))))
;;     (= 1 (in-region-G4-c z a prec))))

;; (define* (in-region-G5? z a #:key (prec 53))
;;   "Check if complex number z is in region G5"
;;   (let ((in-region-G5-c (foreign-library-function libpartition "in_region_G5"
;;                             #:return-type unsigned-short
;;                             #:arg-types (list complex-double double long))))
;;     (= 1 (in-region-G5-c z a prec))))

;; (define* (in-region-G6? z a #:key (prec 53))
;;   "Check if complex number z is in region G6"
;;   (let ((in-region-G6-c (foreign-library-function libpartition "in_region_G6"
;;                             #:return-type unsigned-short
;;                             #:arg-types (list complex-double double long))))
;;     (= 1 (in-region-G6-c z a prec))))

