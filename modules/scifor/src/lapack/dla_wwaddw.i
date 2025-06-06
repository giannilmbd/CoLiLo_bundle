# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/dla_wwaddw.f"
      SUBROUTINE DLA_WWADDW( N, X, Y, W )
*
*     -- LAPACK routine (version 3.2.2)                                 --
*     -- Contributed by James Demmel, Deaglan Halligan, Yozo Hida and --
*     -- Jason Riedy of Univ. of California Berkeley.                 --
*     -- June 2010                                                    --
*
*     -- LAPACK is a software package provided by Univ. of Tennessee, --
*     -- Univ. of California Berkeley and NAG Ltd.                    --
*
      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
      INTEGER            N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * ), W( * )
*     ..
*
*     Purpose
*     =======
*
*     DLA_WWADDW adds a vector W into a doubled-single vector (X, Y).
*
*     This works for all extant IBM's hex and binary floating point
*     arithmetics, but not for decimal.
*
*     Arguments
*     =========
*
*     N      (input) INTEGER
*            The length of vectors X, Y, and W.
*
*     X      (input/output) DOUBLE PRECISION array, dimension (N)
*            The first part of the doubled-single accumulation vector.
*
*     Y      (input/output) DOUBLE PRECISION array, dimension (N)
*            The second part of the doubled-single accumulation vector.
*
*     W      (input) DOUBLE PRECISION array, dimension (N)
*            The vector to be added.
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION   S
      INTEGER            I
*     ..
*     .. Executable Statements ..
*
      DO 10 I = 1, N
        S = X(I) + W(I)
        S = (S + S) - S
        Y(I) = ((X(I) - S) + W(I)) + Y(I)
        X(I) = S
 10   CONTINUE
      RETURN
      END
