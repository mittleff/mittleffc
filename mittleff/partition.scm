(define-module (mittleff partition)
  #:use-module (system foreign)
  #:use-module (system foreign-library)
  #:use-module (mittleff constants)
  #:use-module (srfi srfi-1)
  #:export (compute-r1
            in-region-G0?
            in-region-G1?
            in-region-G2?
            in-region-G3?
            in-region-G4?
            in-region-G5?
            in-region-G6?))

(define libpartition (load-foreign-library
                  "build/libpartition"
                  #:global? #t))

(define* (in-region-G1? z a #:key (acc *default-precision*) (radius *taylor-radius*))
  (let ((c-in-region-G1?
         (pointer->procedure
          int (dynamic-func "in_region_G1" libpartition) `(,complex-double ,double))))
    (if (= 1 (c-in-region-G1? z a)) #t #f)))

(define* (in-region-G2? z a #:key (acc *default-precision*) (radius *taylor-radius*))
  (let ((c-in-region-G2?
         (pointer->procedure
          int (dynamic-func "in_region_G2" libpartition) `(,complex-double ,double))))
    (if (= 1 (c-in-region-G2? z a)) #t #f)))

(define* (in-region-G3? z a #:key (acc *default-precision*) (radius *taylor-radius*))
  (let ((c-in-region-G3?
         (pointer->procedure
          int (dynamic-func "in_region_G3" libpartition) `(,complex-double ,double))))
    (if (= 1 (c-in-region-G3? z a)) #t #f)))

(define* (in-region-G4? z a #:key (acc *default-precision*) (radius *taylor-radius*))
  (let ((c-in-region-G4?
         (pointer->procedure
          int (dynamic-func "in_region_G4" libpartition) `(,complex-double ,double))))
    (if (= 1 (c-in-region-G4? z a)) #t #f)))

(define* (in-region-G5? z a #:key (acc *default-precision*) (radius *taylor-radius*))
  (let ((c-in-region-G5?
         (pointer->procedure
          int (dynamic-func "in_region_G5" libpartition) `(,complex-double ,double ,double))))
    (if (= 1 (c-in-region-G5? z a radius)) #t #f)))

(define* (in-region-G6? z a #:key (acc *default-precision*) (radius *taylor-radius*))
  (let ((c-in-region-G6?
         (pointer->procedure
          int (dynamic-func "in_region_G6" libpartition) `(,complex-double ,double ,double))))
    (if (= 1 (c-in-region-G6? z a radius)) #t #f)))

(define* (compute-r1 a #:key (acc *default-precision*))
  (let ((c-compute-r1
         (pointer->procedure
          double (dynamic-func "compute_r1" libpartition) `(,double ,double))))
    (c-compute-r1 a acc)))
