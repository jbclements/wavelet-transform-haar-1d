#lang typed/racket

;; this file implements a 1-d linear Haar wavelet transform.
;; Tests live in another file, that compares the result of
;; this with a non-fast transform

(require math/array
         typed/rackunit
         plot
         racket/block)

(provide fast-haar-transform)

;; the basic idea here is that--for the Haar basis, at least--all you
;; care about is the sum of some subset of values. This means that after
;; you've extracted all of the high-frequency terms, you can compress
;; by a factor of two without affecting your ability to compute
;; the lower-frequency terms.

;; given a one-dimensional array of even length l and a multiplier k,
;; construct a new one-dimensional array of length l/2 where each term
;; is formed my multiplying k by the sum of a pointwise multiple of
;; pair in the array by the Haar Wavelet pair 1,-1
(: high-freq ((Array Float) Float -> (Array Float)))
(define (high-freq arr k)
  (define old-len (vector-ref (array-shape arr) 0))
  (: new-len Natural)
  (define new-len
    (match (/ old-len 2)
      [(? exact-integer? n) n]
      [else (raise-argument-error 'high-freq "vector of size divisible by 2"
                                  0 arr k)]))
  (for/array: #:shape (vector new-len)
    ([i (in-range new-len)])  : Float
    (ann (* k (- (array-ref arr (vector (* 2 i)))
                 (array-ref arr (vector (add1 (* 2 i))))))
         Float)))

;; given a one-dimensional array of even length l, construct a new
;; vector of length l/2 where each term is the sum of two consecutive
;; terms of the old array
(: squeeze ((Array Float) -> (Array Float)))
(define (squeeze arr)
  (define new-len (/ (vector-ref (array-shape arr) 0) 2))
  (for/array ([i (in-range new-len)]) : Float
    (+ (array-ref arr (vector (* 2 i)))
       (array-ref arr (vector (add1 (* 2 i)))))))

;; given a list of 1-d arrays, produce a new array containing
;; the elements from the arrays in order
(: arrays-join (All (T) ((Listof (Array T)) -> (Array T))))
(define (arrays-join as)
  (define elements
    (foldr (ann sequence-append
                ((Sequenceof T) (Sequenceof T) -> (Sequenceof T)))
           '()
           (map (ann in-array ((Array T) -> (Sequenceof T)))
                as)))
  (for/array ([v elements]) : T
    v))

(define-predicate nonnegative-Float? Nonnegative-Float)
(define-predicate positive-float? Positive-Float)

;; given a 1-d array whose length is a power of two, 
;; use high-freq and squeeze to build a linear-time transform:
(: fast-haar-transform ((Array Float) -> (Array Float)))
(define (fast-haar-transform arr)
  (: orig-len Positive-Float)
  (define orig-len
    (match (exact->inexact (vector-ref (array-shape arr) 0))
      [(? positive-float? n) n]
      [else (error "internal error, length of array not positive-Float")]))
  (define sub-arrays
    (let loop : (Listof (Array Float)) ([arr arr])
      (: len Positive-Float)
      (define len
        (match (exact->inexact (vector-ref (array-shape arr) 0))
          [(? positive-float? n) n]
          [else (error "internal error, length of array not positive-Float")]))
      (: k Float)
      (define k (sqrt
                 (match (/ 1.0 (* orig-len (/ 2.0 len)))
                   [(? positive-float? n) n]
                   [else (error "internal error, expected positive-Float")])))
      (define hf (high-freq arr k))
      (cons hf
            (cond [(= len 2.0)
                   ;; must special-case vector zero
                   (list
                    (list->array
                     (list
                      (* k
                         (+ (array-ref arr (vector 0))
                            (array-ref arr (vector 1)))))))]
                  [else
                   (loop (squeeze arr))]))))
  (arrays-join (reverse sub-arrays)))
