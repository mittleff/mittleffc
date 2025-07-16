(define-module (tests test-utils)
  #:use-module (mittleff utils)
  #:use-module (srfi srfi-1)
  #:use-module (srfi srfi-64)
  #:use-module (ice-9 format)
  #:use-module (ice-9 ftw)
  #:use-module (tests test-approx)
  #:use-module (tests runner))

(test-runner-factory mittleff:test-runner)

(test-begin "accuracy-from-precision")
(test-approx "(acc-from-prec 49)== 1e-15"
	     1e-15 (accuracy-from-prec 49))
            
(test-end "accuracy-from-precision")
