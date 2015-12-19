#lang scribble/manual

@title{Fast 1-D Haar Wavelet Transform}

@author[(author+email "John Clements" "clements@racket-lang.org")]

@(require (for-label racket
                     ))


@defmodule[wavelet]{

                    This package provides functions to perform
 one-dimensional Haar
wavelet transforms in linear time.

@defproc[(fast-haar-transform [array (Array Real)]) (Array Real)]{
Given a one-dimensional array of Reals whose size is a power of two,
return the one-dimensional array of Reals of the same size that represents
the result of the one-dimensional Haar transform.
 }

@defproc[(fast-haar-inverse-transform [array (Array Real)]) (Array Real)]{
 Given a one-dimensional array of Reals whose size is a power of two,
 return the inverse Haar transform. This function is the inverse of
 @racket[fast-haar-transform].}

}