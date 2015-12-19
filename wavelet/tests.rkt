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

  ;; compute the inverse transform the natural way, by
  ;; summing the products of the terms with the basis
  ;; vectors
  (define (slow-haar-inverse-transform arr)
    (for/fold ([arr (make-array (array-shape arr) 0.0)])
              ([v (in-list basis-vectors)]
               [x (in-array arr)])
      (array+ arr (array* v (array x)))))

  ;; double the length of an array by duplicating each
  ;; element
  (define (expand-array arr)
    (for/array ([i (in-range (* 2 (array-size arr)))])
      (array-ref arr (vector (floor (/ i 2))))))

  (check-equal? (expand-array (array #[0.3 -4]))
                (array #[0.3 0.3 -4 -4]))

  ;; use a subset of the elements in a given array of haar coefficients
  ;; to reconstruct a larger array
  (define (haar-reconstruct-step accum coefficients offset mult)
    (define adjustment-array
      (for/array ([i (in-range (array-size accum))])
        (define coefficient
          (* mult
             (array-ref coefficients
                        (vector (+ offset
                                   (floor (/ i 2)))))))
        (cond [(= 0 (modulo i 2)) coefficient]
              [else (- coefficient)])))
    (array+ accum
            adjustment-array))


  (define ((array-= epsilon) a b)
    (and (equal? (array-shape a) (array-shape b))
         (for/and ([aval (in-array a)]
                   [bval (in-array b)])
           (< (abs (- aval bval)) epsilon))))
  
  (check (array-= 1e-10)
         (haar-reconstruct-step (array #[0.3 0.3 -2 -2])
                                (array #[0.1 0.2 0.3 0.4 0.5 0.6])
                                2
                                3.0)
         (array #[1.2 -0.6 -0.8 -3.2]))

  ;; given a 1-d array of reals
  ;; whose length is a power of 2, perform an inverse
  ;; wavelet transform.
  (define (fast-haar-inverse-transform arr)
    (define init-array (make-array #(1)
                                   (* (array-ref arr (vector 0))
                                      (sqrt (/ 1 (array-size arr))))))
    (let loop ([offset 1]
               [accum init-array])
      (cond [(= offset (array-size arr)) accum]
            [else
             (define expanded (expand-array accum))
             (define multiplier (sqrt (/ (array-size accum) (array-size arr))))
             (loop (+ offset (array-size accum))
                   (haar-reconstruct-step expanded
                                          arr offset
                                          multiplier))])))
  
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
   (test-case
    "fast-transform equal to slow-transform"
    (for ([i (in-range POINTS)])
      (check-= (array-ref transformed2 (vector i))
               (array-ref transformed (vector i))
               1e-10)))

   ;; does the inverse slow transform produce the original vector?
   (test-case
    "slow-inverse-transform matches original"
    (define reverse-transformed (slow-haar-inverse-transform transformed2))
    (for ([i (in-range POINTS)])
      (check-= (array-ref reverse-transformed (vector i))
               (array-ref arr (vector i))
               1e-10)))

   (test-case
    "fast-inverse-transform"
    (define reverse-transformed2 (fast-haar-inverse-transform transformed2))
    (for ([i (in-range POINTS)])
      (check-= (array-ref reverse-transformed2 (vector i))
               (array-ref arr (vector i))
               1e-10)))
   ))

