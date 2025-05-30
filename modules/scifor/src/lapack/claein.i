# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/claein.f"
      SUBROUTINE CLAEIN( RIGHTV, NOINIT, N, H, LDH, W, V, B, LDB, RWORK,
     $                   EPS3, SMLNUM, INFO )
*
*  -- LAPACK auxiliary routine (version 3.3.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- April 2011                                                      --
*
*     .. Scalar Arguments ..
      LOGICAL            NOINIT, RIGHTV
      INTEGER            INFO, LDB, LDH, N
      REAL               EPS3, SMLNUM
      COMPLEX            W
*     ..
*     .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            B( LDB, * ), H( LDH, * ), V( * )
*     ..
*
*  Purpose
*  =======
*
*  CLAEIN uses inverse iteration to find a right or left eigenvector
*  corresponding to the eigenvalue W of a complex upper Hessenberg
*  matrix H.
*
*  Arguments
*  =========
*
*  RIGHTV   (input) LOGICAL
*          = .TRUE. : compute right eigenvector;
*          = .FALSE.: compute left eigenvector.
*
*  NOINIT   (input) LOGICAL
*          = .TRUE. : no initial vector supplied in V
*          = .FALSE.: initial vector supplied in V.
*
*  N       (input) INTEGER
*          The order of the matrix H.  N >= 0.
*
*  H       (input) COMPLEX array, dimension (LDH,N)
*          The upper Hessenberg matrix H.
*
*  LDH     (input) INTEGER
*          The leading dimension of the array H.  LDH >= max(1,N).
*
*  W       (input) COMPLEX
*          The eigenvalue of H whose corresponding right or left
*          eigenvector is to be computed.
*
*  V       (input/output) COMPLEX array, dimension (N)
*          On entry, if NOINIT = .FALSE., V must contain a starting
*          vector for inverse iteration; otherwise V need not be set.
*          On exit, V contains the computed eigenvector, normalized so
*          that the component of largest magnitude has magnitude 1; here
*          the magnitude of a complex number (x,y) is taken to be
*          |x| + |y|.
*
*  B       (workspace) COMPLEX array, dimension (LDB,N)
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  RWORK   (workspace) REAL array, dimension (N)
*
*  EPS3    (input) REAL
*          A small machine-dependent value which is used to perturb
*          close eigenvalues, and to replace zero pivots.
*
*  SMLNUM  (input) REAL
*          A machine-dependent value close to the underflow threshold.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          = 1:  inverse iteration did not converge; V is set to the
*                last iterate.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, TENTH
      PARAMETER          ( ONE = 1.0E+0, TENTH = 1.0E-1 )
      COMPLEX            ZERO
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      CHARACTER          NORMIN, TRANS
      INTEGER            I, IERR, ITS, J
      REAL               GROWTO, NRMSML, ROOTN, RTEMP, SCALE, VNORM
      COMPLEX            CDUM, EI, EJ, TEMP, X
*     ..
*     .. External Functions ..
      INTEGER            ICAMAX
      REAL               SCASUM, SCNRM2
      COMPLEX            CLADIV
      EXTERNAL           ICAMAX, SCASUM, SCNRM2, CLADIV
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLATRS, CSSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, MAX, REAL, SQRT
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     GROWTO is the threshold used in the acceptance test for an
*     eigenvector.
*
      ROOTN = SQRT( REAL( N ) )
      GROWTO = TENTH / ROOTN
      NRMSML = MAX( ONE, EPS3*ROOTN )*SMLNUM
*
*     Form B = H - W*I (except that the subdiagonal elements are not
*     stored).
*
      DO 20 J = 1, N
         DO 10 I = 1, J - 1
            B( I, J ) = H( I, J )
   10    CONTINUE
         B( J, J ) = H( J, J ) - W
   20 CONTINUE
*
      IF( NOINIT ) THEN
*
*        Initialize V.
*
         DO 30 I = 1, N
            V( I ) = EPS3
   30    CONTINUE
      ELSE
*
*        Scale supplied initial vector.
*
         VNORM = SCNRM2( N, V, 1 )
         CALL CSSCAL( N, ( EPS3*ROOTN ) / MAX( VNORM, NRMSML ), V, 1 )
      END IF
*
      IF( RIGHTV ) THEN
*
*        LU decomposition with partial pivoting of B, replacing zero
*        pivots by EPS3.
*
         DO 60 I = 1, N - 1
            EI = H( I+1, I )
            IF( CABS1( B( I, I ) ).LT.CABS1( EI ) ) THEN
*
*              Interchange rows and eliminate.
*
               X = CLADIV( B( I, I ), EI )
               B( I, I ) = EI
               DO 40 J = I + 1, N
                  TEMP = B( I+1, J )
                  B( I+1, J ) = B( I, J ) - X*TEMP
                  B( I, J ) = TEMP
   40          CONTINUE
            ELSE
*
*              Eliminate without interchange.
*
               IF( B( I, I ).EQ.ZERO )
     $            B( I, I ) = EPS3
               X = CLADIV( EI, B( I, I ) )
               IF( X.NE.ZERO ) THEN
                  DO 50 J = I + 1, N
                     B( I+1, J ) = B( I+1, J ) - X*B( I, J )
   50             CONTINUE
               END IF
            END IF
   60    CONTINUE
         IF( B( N, N ).EQ.ZERO )
     $      B( N, N ) = EPS3
*
         TRANS = 'N'
*
      ELSE
*
*        UL decomposition with partial pivoting of B, replacing zero
*        pivots by EPS3.
*
         DO 90 J = N, 2, -1
            EJ = H( J, J-1 )
            IF( CABS1( B( J, J ) ).LT.CABS1( EJ ) ) THEN
*
*              Interchange columns and eliminate.
*
               X = CLADIV( B( J, J ), EJ )
               B( J, J ) = EJ
               DO 70 I = 1, J - 1
                  TEMP = B( I, J-1 )
                  B( I, J-1 ) = B( I, J ) - X*TEMP
                  B( I, J ) = TEMP
   70          CONTINUE
            ELSE
*
*              Eliminate without interchange.
*
               IF( B( J, J ).EQ.ZERO )
     $            B( J, J ) = EPS3
               X = CLADIV( EJ, B( J, J ) )
               IF( X.NE.ZERO ) THEN
                  DO 80 I = 1, J - 1
                     B( I, J-1 ) = B( I, J-1 ) - X*B( I, J )
   80             CONTINUE
               END IF
            END IF
   90    CONTINUE
         IF( B( 1, 1 ).EQ.ZERO )
     $      B( 1, 1 ) = EPS3
*
         TRANS = 'C'
*
      END IF
*
      NORMIN = 'N'
      DO 110 ITS = 1, N
*
*        Solve U*x = scale*v for a right eigenvector
*          or U**H *x = scale*v for a left eigenvector,
*        overwriting x on v.
*
         CALL CLATRS( 'Upper', TRANS, 'Nonunit', NORMIN, N, B, LDB, V,
     $                SCALE, RWORK, IERR )
         NORMIN = 'Y'
*
*        Test for sufficient growth in the norm of v.
*
         VNORM = SCASUM( N, V, 1 )
         IF( VNORM.GE.GROWTO*SCALE )
     $      GO TO 120
*
*        Choose new orthogonal starting vector and try again.
*
         RTEMP = EPS3 / ( ROOTN+ONE )
         V( 1 ) = EPS3
         DO 100 I = 2, N
            V( I ) = RTEMP
  100    CONTINUE
         V( N-ITS+1 ) = V( N-ITS+1 ) - EPS3*ROOTN
  110 CONTINUE
*
*     Failure to find eigenvector in N iterations.
*
      INFO = 1
*
  120 CONTINUE
*
*     Normalize eigenvector.
*
      I = ICAMAX( N, V, 1 )
      CALL CSSCAL( N, ONE / CABS1( V( I ) ), V, 1 )
*
      RETURN
*
*     End of CLAEIN
*
      END
