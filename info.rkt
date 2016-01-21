#lang setup/infotab

(define collection 'multi)

(define deps
  (list
   "base"
   "math-lib"
   "plot-gui-lib"
   "typed-racket-lib"
   "typed-racket-more"))

(define build-deps
  (list "rackunit-lib"
        "racket-doc"
        "scribble-lib"))

