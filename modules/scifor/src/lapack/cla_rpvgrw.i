# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/cla_rpvgrw.f"
      REAL FUNCTION CLA_RPVGRW( N, NCOLS, A, LDA, AF, LDAF )
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
      INTEGER            N, NCOLS, LDA, LDAF
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), AF( LDAF, * )
*     ..
*
*  Purpose
*  =======
*
*  CLA_RPVGRW computes the reciprocal pivot growth factor
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
*     NCOLS   (input) INTEGER
*     The number of columns of the matrix A. NCOLS >= 0.
*
*     A       (input) COMPLEX array, dimension (LDA,N)
*     On entry, the N-by-N matrix A.
*
*     LDA     (input) INTEGER
*     The leading dimension of the array A.  LDA >= max(1,N).
*
*     AF      (input) COMPLEX array, dimension (LDAF,N)
*     The factors L and U from the factorization
*     A = P*L*U as computed by CGETRF.
*
*     LDAF    (input) INTEGER
*     The leading dimension of the array AF.  LDAF >= max(1,N).
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
      REAL               AMAX, UMAX, RPVGRW
      COMPLEX            ZDUM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, ABS, REAL, AIMAG
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
      RPVGRW = 1.0

      DO J = 1, NCOLS
         AMAX = 0.0
         UMAX = 0.0
         DO I = 1, N
            AMAX = MAX( CABS1( A( I, J ) ), AMAX )
         END DO
         DO I = 1, J
            UMAX = MAX( CABS1( AF( I, J ) ), UMAX )
         END DO
         IF ( UMAX /= 0.0 ) THEN
            RPVGRW = MIN( AMAX / UMAX, RPVGRW )
         END IF
      END DO
      CLA_RPVGRW = RPVGRW
      END
