# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/cptrfs.f"
      SUBROUTINE CPTRFS( UPLO, N, NRHS, D, E, DF, EF, B, LDB, X, LDX,
     $                   FERR, BERR, WORK, RWORK, INFO )
*
*  -- LAPACK routine (version 3.3.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- April 2011                                                      --
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDB, LDX, N, NRHS
*     ..
*     .. Array Arguments ..
      REAL               BERR( * ), D( * ), DF( * ), FERR( * ),
     $                   RWORK( * )
      COMPLEX            B( LDB, * ), E( * ), EF( * ), WORK( * ),
     $                   X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  CPTRFS improves the computed solution to a system of linear
*  equations when the coefficient matrix is Hermitian positive definite
*  and tridiagonal, and provides error bounds and backward error
*  estimates for the solution.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the superdiagonal or the subdiagonal of the
*          tridiagonal matrix A is stored and the form of the
*          factorization:
*          = 'U':  E is the superdiagonal of A, and A = U**H*D*U;
*          = 'L':  E is the subdiagonal of A, and A = L*D*L**H.
*          (The two forms are equivalent if A is real.)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  D       (input) REAL array, dimension (N)
*          The n real diagonal elements of the tridiagonal matrix A.
*
*  E       (input) COMPLEX array, dimension (N-1)
*          The (n-1) off-diagonal elements of the tridiagonal matrix A
*          (see UPLO).
*
*  DF      (input) REAL array, dimension (N)
*          The n diagonal elements of the diagonal matrix D from
*          the factorization computed by CPTTRF.
*
*  EF      (input) COMPLEX array, dimension (N-1)
*          The (n-1) off-diagonal elements of the unit bidiagonal
*          factor U or L from the factorization computed by CPTTRF
*          (see UPLO).
*
*  B       (input) COMPLEX array, dimension (LDB,NRHS)
*          The right hand side matrix B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  X       (input/output) COMPLEX array, dimension (LDX,NRHS)
*          On entry, the solution matrix X, as computed by CPTTRS.
*          On exit, the improved solution matrix X.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.  LDX >= max(1,N).
*
*  FERR    (output) REAL array, dimension (NRHS)
*          The forward error bound for each solution vector
*          X(j) (the j-th column of the solution matrix X).
*          If XTRUE is the true solution corresponding to X(j), FERR(j)
*          is an estimated upper bound for the magnitude of the largest
*          element in (X(j) - XTRUE) divided by the magnitude of the
*          largest element in X(j).
*
*  BERR    (output) REAL array, dimension (NRHS)
*          The componentwise relative backward error of each solution
*          vector X(j) (i.e., the smallest relative change in
*          any element of A or B that makes X(j) an exact solution).
*
*  WORK    (workspace) COMPLEX array, dimension (N)
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
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0E+0 )
      REAL               THREE
      PARAMETER          ( THREE = 3.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            COUNT, I, IX, J, NZ
      REAL               EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN
      COMPLEX            BI, CX, DX, EX, ZDUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ISAMAX
      REAL               SLAMCH
      EXTERNAL           LSAME, ISAMAX, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CAXPY, CPTTRS, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CMPLX, CONJG, MAX, REAL
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
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CPTRFS', -INFO )
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
      NZ = 4
      EPS = SLAMCH( 'Epsilon' )
      SAFMIN = SLAMCH( 'Safe minimum' )
      SAFE1 = NZ*SAFMIN
      SAFE2 = SAFE1 / EPS
*
*     Do for each right hand side
*
      DO 100 J = 1, NRHS
*
         COUNT = 1
         LSTRES = THREE
   20    CONTINUE
*
*        Loop until stopping criterion is satisfied.
*
*        Compute residual R = B - A * X.  Also compute
*        abs(A)*abs(x) + abs(b) for use in the backward error bound.
*
         IF( UPPER ) THEN
            IF( N.EQ.1 ) THEN
               BI = B( 1, J )
               DX = D( 1 )*X( 1, J )
               WORK( 1 ) = BI - DX
               RWORK( 1 ) = CABS1( BI ) + CABS1( DX )
            ELSE
               BI = B( 1, J )
               DX = D( 1 )*X( 1, J )
               EX = E( 1 )*X( 2, J )
               WORK( 1 ) = BI - DX - EX
               RWORK( 1 ) = CABS1( BI ) + CABS1( DX ) +
     $                      CABS1( E( 1 ) )*CABS1( X( 2, J ) )
               DO 30 I = 2, N - 1
                  BI = B( I, J )
                  CX = CONJG( E( I-1 ) )*X( I-1, J )
                  DX = D( I )*X( I, J )
                  EX = E( I )*X( I+1, J )
                  WORK( I ) = BI - CX - DX - EX
                  RWORK( I ) = CABS1( BI ) +
     $                         CABS1( E( I-1 ) )*CABS1( X( I-1, J ) ) +
     $                         CABS1( DX ) + CABS1( E( I ) )*
     $                         CABS1( X( I+1, J ) )
   30          CONTINUE
               BI = B( N, J )
               CX = CONJG( E( N-1 ) )*X( N-1, J )
               DX = D( N )*X( N, J )
               WORK( N ) = BI - CX - DX
               RWORK( N ) = CABS1( BI ) + CABS1( E( N-1 ) )*
     $                      CABS1( X( N-1, J ) ) + CABS1( DX )
            END IF
         ELSE
            IF( N.EQ.1 ) THEN
               BI = B( 1, J )
               DX = D( 1 )*X( 1, J )
               WORK( 1 ) = BI - DX
               RWORK( 1 ) = CABS1( BI ) + CABS1( DX )
            ELSE
               BI = B( 1, J )
               DX = D( 1 )*X( 1, J )
               EX = CONJG( E( 1 ) )*X( 2, J )
               WORK( 1 ) = BI - DX - EX
               RWORK( 1 ) = CABS1( BI ) + CABS1( DX ) +
     $                      CABS1( E( 1 ) )*CABS1( X( 2, J ) )
               DO 40 I = 2, N - 1
                  BI = B( I, J )
                  CX = E( I-1 )*X( I-1, J )
                  DX = D( I )*X( I, J )
                  EX = CONJG( E( I ) )*X( I+1, J )
                  WORK( I ) = BI - CX - DX - EX
                  RWORK( I ) = CABS1( BI ) +
     $                         CABS1( E( I-1 ) )*CABS1( X( I-1, J ) ) +
     $                         CABS1( DX ) + CABS1( E( I ) )*
     $                         CABS1( X( I+1, J ) )
   40          CONTINUE
               BI = B( N, J )
               CX = E( N-1 )*X( N-1, J )
               DX = D( N )*X( N, J )
               WORK( N ) = BI - CX - DX
               RWORK( N ) = CABS1( BI ) + CABS1( E( N-1 ) )*
     $                      CABS1( X( N-1, J ) ) + CABS1( DX )
            END IF
         END IF
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
         S = ZERO
         DO 50 I = 1, N
            IF( RWORK( I ).GT.SAFE2 ) THEN
               S = MAX( S, CABS1( WORK( I ) ) / RWORK( I ) )
            ELSE
               S = MAX( S, ( CABS1( WORK( I ) )+SAFE1 ) /
     $             ( RWORK( I )+SAFE1 ) )
            END IF
   50    CONTINUE
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
            CALL CPTTRS( UPLO, N, 1, DF, EF, WORK, N, INFO )
            CALL CAXPY( N, CMPLX( ONE ), WORK, 1, X( 1, J ), 1 )
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
         DO 60 I = 1, N
            IF( RWORK( I ).GT.SAFE2 ) THEN
               RWORK( I ) = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I )
            ELSE
               RWORK( I ) = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I ) +
     $                      SAFE1
            END IF
   60    CONTINUE
         IX = ISAMAX( N, RWORK, 1 )
         FERR( J ) = RWORK( IX )
*
*        Estimate the norm of inv(A).
*
*        Solve M(A) * x = e, where M(A) = (m(i,j)) is given by
*
*           m(i,j) =  abs(A(i,j)), i = j,
*           m(i,j) = -abs(A(i,j)), i .ne. j,
*
*        and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**H.
*
*        Solve M(L) * x = e.
*
         RWORK( 1 ) = ONE
         DO 70 I = 2, N
            RWORK( I ) = ONE + RWORK( I-1 )*ABS( EF( I-1 ) )
   70    CONTINUE
*
*        Solve D * M(L)**H * x = b.
*
         RWORK( N ) = RWORK( N ) / DF( N )
         DO 80 I = N - 1, 1, -1
            RWORK( I ) = RWORK( I ) / DF( I ) +
     $                   RWORK( I+1 )*ABS( EF( I ) )
   80    CONTINUE
*
*        Compute norm(inv(A)) = max(x(i)), 1<=i<=n.
*
         IX = ISAMAX( N, RWORK, 1 )
         FERR( J ) = FERR( J )*ABS( RWORK( IX ) )
*
*        Normalize error.
*
         LSTRES = ZERO
         DO 90 I = 1, N
            LSTRES = MAX( LSTRES, ABS( X( I, J ) ) )
   90    CONTINUE
         IF( LSTRES.NE.ZERO )
     $      FERR( J ) = FERR( J ) / LSTRES
*
  100 CONTINUE
*
      RETURN
*
*     End of CPTRFS
*
      END
