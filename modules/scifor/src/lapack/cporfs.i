# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/cporfs.f"
      SUBROUTINE CPORFS( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X,
     $                   LDX, FERR, BERR, WORK, RWORK, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     Modified to call CLACN2 in place of CLACON, 10 Feb 03, SJH.
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
*     ..
*     .. Array Arguments ..
      REAL               BERR( * ), FERR( * ), RWORK( * )
      COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
     $                   WORK( * ), X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  CPORFS improves the computed solution to a system of linear
*  equations when the coefficient matrix is Hermitian positive definite,
*  and provides error bounds and backward error estimates for the
*  solution.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrices B and X.  NRHS >= 0.
*
*  A       (input) COMPLEX array, dimension (LDA,N)
*          The Hermitian matrix A.  If UPLO = 'U', the leading N-by-N
*          upper triangular part of A contains the upper triangular part
*          of the matrix A, and the strictly lower triangular part of A
*          is not referenced.  If UPLO = 'L', the leading N-by-N lower
*          triangular part of A contains the lower triangular part of
*          the matrix A, and the strictly upper triangular part of A is
*          not referenced.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  AF      (input) COMPLEX array, dimension (LDAF,N)
*          The triangular factor U or L from the Cholesky factorization
*          A = U**H*U or A = L*L**H, as computed by CPOTRF.
*
*  LDAF    (input) INTEGER
*          The leading dimension of the array AF.  LDAF >= max(1,N).
*
*  B       (input) COMPLEX array, dimension (LDB,NRHS)
*          The right hand side matrix B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  X       (input/output) COMPLEX array, dimension (LDX,NRHS)
*          On entry, the solution matrix X, as computed by CPOTRS.
*          On exit, the improved solution matrix X.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.  LDX >= max(1,N).
*
*  FERR    (output) REAL array, dimension (NRHS)
*          The estimated forward error bound for each solution vector
*          X(j) (the j-th column of the solution matrix X).
*          If XTRUE is the true solution corresponding to X(j), FERR(j)
*          is an estimated upper bound for the magnitude of the largest
*          element in (X(j) - XTRUE) divided by the magnitude of the
*          largest element in X(j).  The estimate is as reliable as
*          the estimate for RCOND, and is almost always a slight
*          overestimate of the true error.
*
*  BERR    (output) REAL array, dimension (NRHS)
*          The componentwise relative backward error of each solution
*          vector X(j) (i.e., the smallest relative change in
*          any element of A or B that makes X(j) an exact solution).
*
*  WORK    (workspace) COMPLEX array, dimension (2*N)
*
*  RWORK   (workspace) REAL array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Internal Parameters
*  ===================
*
*  ITMAX is the maximum number of steps of iterative refinement.
*
*  ====================================================================
*
*     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
      REAL               TWO
      PARAMETER          ( TWO = 2.0E+0 )
      REAL               THREE
      PARAMETER          ( THREE = 3.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            COUNT, I, J, K, KASE, NZ
      REAL               EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK
      COMPLEX            ZDUM
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           CAXPY, CCOPY, CHEMV, CLACN2, CPOTRS, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, MAX, REAL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH
      EXTERNAL           LSAME, SLAMCH
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDAF.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CPORFS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) THEN
         DO 10 J = 1, NRHS
            FERR( J ) = ZERO
            BERR( J ) = ZERO
   10    CONTINUE
         RETURN
      END IF
*
*     NZ = maximum number of nonzero elements in each row of A, plus 1
*
      NZ = N + 1
      EPS = SLAMCH( 'Epsilon' )
      SAFMIN = SLAMCH( 'Safe minimum' )
      SAFE1 = NZ*SAFMIN
      SAFE2 = SAFE1 / EPS
*
*     Do for each right hand side
*
      DO 140 J = 1, NRHS
*
         COUNT = 1
         LSTRES = THREE
   20    CONTINUE
*
*        Loop until stopping criterion is satisfied.
*
*        Compute residual R = B - A * X
*
         CALL CCOPY( N, B( 1, J ), 1, WORK, 1 )
         CALL CHEMV( UPLO, N, -ONE, A, LDA, X( 1, J ), 1, ONE, WORK, 1 )
*
*        Compute componentwise relative backward error from formula
*
*        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )
*
*        where abs(Z) is the componentwise absolute value of the matrix
*        or vector Z.  If the i-th component of the denominator is less
*        than SAFE2, then SAFE1 is added to the i-th components of the
*        numerator and denominator before dividing.
*
         DO 30 I = 1, N
            RWORK( I ) = CABS1( B( I, J ) )
   30    CONTINUE
*
*        Compute abs(A)*abs(X) + abs(B).
*
         IF( UPPER ) THEN
            DO 50 K = 1, N
               S = ZERO
               XK = CABS1( X( K, J ) )
               DO 40 I = 1, K - 1
                  RWORK( I ) = RWORK( I ) + CABS1( A( I, K ) )*XK
                  S = S + CABS1( A( I, K ) )*CABS1( X( I, J ) )
   40          CONTINUE
               RWORK( K ) = RWORK( K ) + ABS( REAL( A( K, K ) ) )*XK + S
   50       CONTINUE
         ELSE
            DO 70 K = 1, N
               S = ZERO
               XK = CABS1( X( K, J ) )
               RWORK( K ) = RWORK( K ) + ABS( REAL( A( K, K ) ) )*XK
               DO 60 I = K + 1, N
                  RWORK( I ) = RWORK( I ) + CABS1( A( I, K ) )*XK
                  S = S + CABS1( A( I, K ) )*CABS1( X( I, J ) )
   60          CONTINUE
               RWORK( K ) = RWORK( K ) + S
   70       CONTINUE
         END IF
         S = ZERO
         DO 80 I = 1, N
            IF( RWORK( I ).GT.SAFE2 ) THEN
               S = MAX( S, CABS1( WORK( I ) ) / RWORK( I ) )
            ELSE
               S = MAX( S, ( CABS1( WORK( I ) )+SAFE1 ) /
     $             ( RWORK( I )+SAFE1 ) )
            END IF
   80    CONTINUE
         BERR( J ) = S
*
*        Test stopping criterion. Continue iterating if
*           1) The residual BERR(J) is larger than machine epsilon, and
*           2) BERR(J) decreased by at least a factor of 2 during the
*              last iteration, and
*           3) At most ITMAX iterations tried.
*
         IF( BERR( J ).GT.EPS .AND. TWO*BERR( J ).LE.LSTRES .AND.
     $       COUNT.LE.ITMAX ) THEN
*
*           Update solution and try again.
*
            CALL CPOTRS( UPLO, N, 1, AF, LDAF, WORK, N, INFO )
            CALL CAXPY( N, ONE, WORK, 1, X( 1, J ), 1 )
            LSTRES = BERR( J )
            COUNT = COUNT + 1
            GO TO 20
         END IF
*
*        Bound error from formula
*
*        norm(X - XTRUE) / norm(X) .le. FERR =
*        norm( abs(inv(A))*
*           ( abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) / norm(X)
*
*        where
*          norm(Z) is the magnitude of the largest component of Z
*          inv(A) is the inverse of A
*          abs(Z) is the componentwise absolute value of the matrix or
*             vector Z
*          NZ is the maximum number of nonzeros in any row of A, plus 1
*          EPS is machine epsilon
*
*        The i-th component of abs(R)+NZ*EPS*(abs(A)*abs(X)+abs(B))
*        is incremented by SAFE1 if the i-th component of
*        abs(A)*abs(X) + abs(B) is less than SAFE2.
*
*        Use CLACN2 to estimate the infinity-norm of the matrix
*           inv(A) * diag(W),
*        where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) )))
*
         DO 90 I = 1, N
            IF( RWORK( I ).GT.SAFE2 ) THEN
               RWORK( I ) = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I )
            ELSE
               RWORK( I ) = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I ) +
     $                      SAFE1
            END IF
   90    CONTINUE
*
         KASE = 0
  100    CONTINUE
         CALL CLACN2( N, WORK( N+1 ), WORK, FERR( J ), KASE, ISAVE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.1 ) THEN
*
*              Multiply by diag(W)*inv(A**H).
*
               CALL CPOTRS( UPLO, N, 1, AF, LDAF, WORK, N, INFO )
               DO 110 I = 1, N
                  WORK( I ) = RWORK( I )*WORK( I )
  110          CONTINUE
            ELSE IF( KASE.EQ.2 ) THEN
*
*              Multiply by inv(A)*diag(W).
*
               DO 120 I = 1, N
                  WORK( I ) = RWORK( I )*WORK( I )
  120          CONTINUE
               CALL CPOTRS( UPLO, N, 1, AF, LDAF, WORK, N, INFO )
            END IF
            GO TO 100
         END IF
*
*        Normalize error.
*
         LSTRES = ZERO
         DO 130 I = 1, N
            LSTRES = MAX( LSTRES, CABS1( X( I, J ) ) )
  130    CONTINUE
         IF( LSTRES.NE.ZERO )
     $      FERR( J ) = FERR( J ) / LSTRES
*
  140 CONTINUE
*
      RETURN
*
*     End of CPORFS
*
      END
