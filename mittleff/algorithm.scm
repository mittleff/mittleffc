(define-module (mittleff algorithm)
  #:use-module (srfi srfi-1)
  #:use-module (system foreign)
  #:use-module (system foreign-library)
  #:use-module (mittleff rgamma)
  #:use-module (mittleff erfc)
  #:use-module (mittleff quad)
  #:use-module (mittleff macros)
  #:use-module (mittleff constants)
  #:use-module (mittleff partition)
  #:export (taylor
            asymptotics
            integral-rep))

(define libalgorithm (load-foreign-library
                      "build/libalgorithm"
                      #:global? #t))

(define* (taylor a b z #:key (acc *default-precision*))
  (let ((series
         (pointer->procedure
          complex-double (dynamic-func "series" libalgorithm) `(,double ,double ,complex-double ,double ,int))))
    (series a b z acc 0)))

(define* (asymptotics a b z #:key (acc *default-precision*))
  (let ((asymptoticss
         (pointer->procedure
          complex-double (dynamic-func "asymptotics" libalgorithm) `(,double ,double ,complex-double ,double))))
    (asymptoticss a b z acc)))

(define fn-a
  (pointer->procedure
   complex-double (dynamic-func "fn_a" libalgorithm) `(,complex-double ,double ,double ,complex-double)))

(define omega (pointer->procedure
               complex-double (dynamic-func "omega" libalgorithm) `(,complex-double ,complex-double ,double ,double)))

(define fn-b
  (pointer->procedure
   complex-double (dynamic-func "fn_b" libalgorithm) `(,double ,double ,double ,complex-double ,double)))

(define fn-c
  (pointer->procedure
   complex-double (dynamic-func "fn_c" libalgorithm) `(,double ,double ,double ,complex-double ,double)))

(define* (integral-rep a b z  #:key (region 5) (acc *default-precision*) (radius *taylor-radius*))
  (integral-rep-aux a b z #:region
                    (cond
                     ((in-region-G5? z a #:acc acc #:radius radius) 5)
                     ((in-region-G6? z a #:acc acc #:radius radius) 6)) #:acc acc))

(define* (integral-rep-aux a b z  #:key (region 5) (acc *default-precision*))
  (define Phi (pointer->procedure
               double (dynamic-func "Phi_i" libalgorithm) `(,int ,double)))
  (define Rho (pointer->procedure
               double (dynamic-func "Rho_beta" libalgorithm) `(,double)))
  (define Rib (pointer->procedure
               double (dynamic-func "Rib" libalgorithm) `(,int ,double ,double ,complex-double ,double)))
  (let ((Ai (cond ((= region 5) (fn-a z a b 0)) ((= region 6) 0)))
        (Bi (quad (lambda (r) (fn-b r a b z (Phi region a))) (Rho b) (Rib region a b z acc)))
        (Ci (if (> b 1) (quad (lambda (ph) (fn-c ph a b z 1/2)) (- (Phi region a)) (Phi region a)) 0)))
    ;;(format #t "~a | ~a | ~a\n" Ai Bi Ci)
    (+ Ai Bi Ci)))

