# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/cla_syamv.f"
      SUBROUTINE CLA_SYAMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y,
     $                      INCY )
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
      REAL               ALPHA, BETA
      INTEGER            INCX, INCY, LDA, N
      INTEGER            UPLO
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), X( * )
      REAL               Y( * )
*     ..
*
*  Purpose
*  =======
*
*  CLA_SYAMV  performs the matrix-vector operation
*
*          y := alpha*abs(A)*abs(x) + beta*abs(y),
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  n by n symmetric matrix.
*
*  This function is primarily used in calculating error bounds.
*  To protect against underflow during evaluation, components in
*  the resulting vector are perturbed away from zero by (N+1)
*  times the underflow threshold.  To prevent unnecessarily large
*  errors for block-structure embedded in general matrices,
*  "symbolically" zero components are not perturbed.  A zero
*  entry is considered "symbolic" if all multiplications involved
*  in computing that entry have at least one zero multiplicand.
*
*  Arguments
*  ==========
*
*  UPLO    (input) INTEGER
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows:
*
*              UPLO = BLAS_UPPER   Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = BLAS_LOWER   Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N       (input) INTEGER
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA   (input) REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - COMPLEX             array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA     (input) INTEGER
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X       (input) COMPLEX array, dimension
*           ( 1 + ( n - 1 )*abs( INCX ) )
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX    (input) INTEGER
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA    (input) REAL            .
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y       (input/output) REAL array, dimension
*           ( 1 + ( n - 1 )*abs( INCY ) )
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY    (input) INTEGER
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  Further Details
*  ===============
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*  -- Modified for the absolute-value product, April 2006
*     Jason Riedy, UC Berkeley
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            SYMB_ZERO
      REAL               TEMP, SAFE1
      INTEGER            I, INFO, IY, J, JX, KX, KY
      COMPLEX            ZDUM
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, SLAMCH
      REAL               SLAMCH
*     ..
*     .. External Functions ..
      EXTERNAL           ILAUPLO
      INTEGER            ILAUPLO
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, ABS, SIGN, REAL, AIMAG
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( REAL ( ZDUM ) ) + ABS( AIMAG ( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( UPLO.NE.ILAUPLO( 'U' ) .AND.
     $         UPLO.NE.ILAUPLO( 'L' ) )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 5
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 7
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SSYMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set up the start points in  X  and  Y.
*
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 )*INCY
      END IF
*
*     Set SAFE1 essentially to be the underflow threshold times the
*     number of additions in each row.
*
      SAFE1 = SLAMCH( 'Safe minimum' )
      SAFE1 = (N+1)*SAFE1
*
*     Form  y := alpha*abs(A)*abs(x) + beta*abs(y).
*
*     The O(N^2) SYMB_ZERO tests could be replaced by O(N) queries to
*     the inexact flag.  Still doesn't help change the iteration order
*     to per-column.
*
      IY = KY
      IF ( INCX.EQ.1 ) THEN
         IF ( UPLO .EQ. ILAUPLO( 'U' ) ) THEN
            DO I = 1, N
               IF ( BETA .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0
               ELSE IF ( Y( IY ) .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
               ELSE
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               END IF
               IF ( ALPHA .NE. ZERO ) THEN
                  DO J = 1, I
                     TEMP = CABS1( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND.
     $                    ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP
                  END DO
                  DO J = I+1, N
                     TEMP = CABS1( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND.
     $                    ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP
                  END DO
               END IF

               IF ( .NOT.SYMB_ZERO )
     $              Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         ELSE
            DO I = 1, N
               IF ( BETA .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0
               ELSE IF ( Y( IY ) .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
               ELSE
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               END IF
               IF ( ALPHA .NE. ZERO ) THEN
                  DO J = 1, I
                     TEMP = CABS1( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND.
     $                    ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP
                  END DO
                  DO J = I+1, N
                     TEMP = CABS1( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND.
     $                    ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP
                  END DO
               END IF

               IF ( .NOT.SYMB_ZERO )
     $              Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         END IF
      ELSE
         IF ( UPLO .EQ. ILAUPLO( 'U' ) ) THEN
            DO I = 1, N
               IF ( BETA .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0
               ELSE IF ( Y( IY ) .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
               ELSE
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               END IF
               JX = KX
               IF ( ALPHA .NE. ZERO ) THEN
                  DO J = 1, I
                     TEMP = CABS1( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND.
     $                    ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
                  DO J = I+1, N
                     TEMP = CABS1( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND.
     $                    ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
               END IF

               IF ( .NOT.SYMB_ZERO )
     $              Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         ELSE
            DO I = 1, N
               IF ( BETA .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0
               ELSE IF ( Y( IY ) .EQ. ZERO ) THEN
                  SYMB_ZERO = .TRUE.
               ELSE
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               END IF
               JX = KX
               IF ( ALPHA .NE. ZERO ) THEN
                  DO J = 1, I
                     TEMP = CABS1( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND.
     $                    ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
                  DO J = I+1, N
                     TEMP = CABS1( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND.
     $                    ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
               END IF

               IF ( .NOT.SYMB_ZERO )
     $              Y( IY ) = Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         END IF

      END IF
*
      RETURN
*
*     End of CLA_SYAMV
*
      END
