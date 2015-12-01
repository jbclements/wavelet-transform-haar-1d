#lang racket

(module+ test
  (require rsound
           plot
           math/array
           rackunit
           racket/block
           "haar-1d.rkt")
  
  
  (define FR 44100)
  (define (s x) (* FR x))
  
  ;; basic dimension for testing is 2^N
  (define N 5)
  (define POINTS (expt 2 N))
  
  ;; build a list of basis vectors of dimension 2^N
  
  (define dc-vector
    (make-array (vector POINTS) (sqrt (/ 1.0 POINTS))))
  
  (define (build-basis-vector i j)
    (define livelen (expt 2 (- N i)))
    (define up (/ livelen 2))
    (define offset (* j livelen))
    (define val (sqrt (/ 1.0 livelen)))
    (build-array
     (vector POINTS)
     (λ (idx)
       (define k (vector-ref idx 0))
       (cond [(< k offset) 0.0]
             [(< k (+ offset up)) val]
             [(< k (+ offset livelen)) (- val)]
             [else 0.0]))))
  
  (define basis-vectors
    (cons dc-vector
          (for*/list ([i (in-range 0 N)]
                      [j (in-range (expt 2 i))])
            (build-basis-vector i j))))
  
  ;; compute the dot product of two vectors. Note no conjugation!
  (define (dot-product a b)
    (array-ref (array-axis-fold (array* a b) 0 +) #()))
  
  ;; ensure that the basis is orthonormal:
  ;; this is a good test, but it takes a long long time for big
  ;; values of N (say, > 7)
  (for* ([a (in-list basis-vectors)]
         [b (in-list basis-vectors)])
    (cond [(eq? a b)
           (check-= (dot-product a b) 1.0 1e-10)]
          [else
           (check-= (dot-product a b) 0.0 1e-10)]))
  
  ;; compute the transform the natural way, by taking
  ;; the dot product with each of the basis vectors:
  (define (slow-haar-transform arr)
    (for/list ([v (in-list basis-vectors)])
      (array-axis-fold (array* v arr) 0 +)))
  
  (block
   (random-seed 779287)
   ;; an array of random reals:
   (define arr
     (build-array (vector POINTS)
                  (λ (idx) (random))))
   ;; transformed using n^2 technique & basis vectors:
   (define transformed
     (for/array ([v (in-list basis-vectors)])
       (array-ref (array-axis-fold (array* v arr) 0 +) #())))
   
   ;; tests of high-freq and squeeze commented out because these
   ;; aren't exported by the typed file:
   #;(;; transformed using high-freq:
      (define hf (high-freq arr (sqrt 1/2)))
      ;; should match top half of terms in transformed
      (for ([i (in-range (/ POINTS 2))])
        (check-= (array-ref hf (vector i))
                 (array-ref transformed (vector (+ (/ POINTS 2) i)))
                 1e-10))
      ;; now squeeze and use high-freq on the result:
      (define squeezed (squeeze arr))
      (define mf (high-freq squeezed (sqrt 1/4)))
      ;; should match second quarter of terms in transformed
      (for ([i (in-range (/ POINTS 4))])
        (check-= (array-ref mf (vector i))
                 (array-ref transformed (vector (+ (/ POINTS 4) i)))
                 1e-10)))
   ;; do the whole thing at once:
   (define transformed2 (fast-haar-transform arr))
   (for ([i (in-range POINTS)])
     (check-= (array-ref transformed (vector i))
              (array-ref transformed2 (vector i))
              1e-10))
   ))


