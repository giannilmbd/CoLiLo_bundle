# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/dla_gbrpvgrw.f"
      DOUBLE PRECISION FUNCTION DLA_GBRPVGRW( N, KL, KU, NCOLS, AB,
     $                                        LDAB, AFB, LDAFB )
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
      INTEGER            N, KL, KU, NCOLS, LDAB, LDAFB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   AB( LDAB, * ), AFB( LDAFB, * )
*     ..
*
*  Purpose
*  =======
*
*  DLA_GBRPVGRW computes the reciprocal pivot growth factor
*  norm(A)/norm(U). The "max absolute element" norm is used. If this is
*  much less than 1, the stability of the LU factorization of the
*  (equilibrated) matrix A could be poor. This also means that the
*  solution X, estimated condition numbers, and error bounds could be
*  unreliable.
*
*  Arguments
*  =========
*
*     N       (input) INTEGER
*     The number of linear equations, i.e., the order of the
*     matrix A.  N >= 0.
*
*     KL      (input) INTEGER
*     The number of subdiagonals within the band of A.  KL >= 0.
*
*     KU      (input) INTEGER
*     The number of superdiagonals within the band of A.  KU >= 0.
*
*     NCOLS   (input) INTEGER
*     The number of columns of the matrix A.  NCOLS >= 0.
*
*     AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
*     On entry, the matrix A in band storage, in rows 1 to KL+KU+1.
*     The j-th column of A is stored in the j-th column of the
*     array AB as follows:
*     AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)
*
*     LDAB    (input) INTEGER
*     The leading dimension of the array AB.  LDAB >= KL+KU+1.
*
*     AFB     (input) DOUBLE PRECISION array, dimension (LDAFB,N)
*     Details of the LU factorization of the band matrix A, as
*     computed by DGBTRF.  U is stored as an upper triangular
*     band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,
*     and the multipliers used during the factorization are stored
*     in rows KL+KU+2 to 2*KL+KU+1.
*
*     LDAFB   (input) INTEGER
*     The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J, KD
      DOUBLE PRECISION   AMAX, UMAX, RPVGRW
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      RPVGRW = 1.0D+0

      KD = KU + 1
      DO J = 1, NCOLS
         AMAX = 0.0D+0
         UMAX = 0.0D+0
         DO I = MAX( J-KU, 1 ), MIN( J+KL, N )
            AMAX = MAX( ABS( AB( KD+I-J, J)), AMAX )
         END DO
         DO I = MAX( J-KU, 1 ), J
            UMAX = MAX( ABS( AFB( KD+I-J, J ) ), UMAX )
         END DO
         IF ( UMAX /= 0.0D+0 ) THEN
            RPVGRW = MIN( AMAX / UMAX, RPVGRW )
         END IF
      END DO
      DLA_GBRPVGRW = RPVGRW
      END
