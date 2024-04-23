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

(define (fn-a z a b x)
  (* (/ a) (expt z (/ (- 1 b) a)) (exp (* (expt z (/ a)) (cos (/ x a))))))

(define (omega x y a b)
  (+ (* (expt x (/ a)) (sin (/ y a))) (* y (+ 1 (/ (- 1 b) a)))))

(define (fn-b r a b z phi)
  (let ((w (omega r phi a b)))
    (* (/ pi) (fn-a r a b phi) (/ (- (* r (sin (- w phi))) (* z (sin w)))
                                  (+ (expt r 2) (* -2 r z (cos phi)) (expt z 2))))))

(define (fn-c ph a b z rho)
  (let ((w (omega rho ph a b)))
    (* (/ rho (* 2 pi)) (fn-a rho a b ph) (/ (exp (* J w)) (- (* rho (exp (* J ph))) z)))))

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
    (+ Ai Bi Ci)))

