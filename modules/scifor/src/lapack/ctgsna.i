# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/ctgsna.f"
      SUBROUTINE CTGSNA( JOB, HOWMNY, SELECT, N, A, LDA, B, LDB, VL,
     $                   LDVL, VR, LDVR, S, DIF, MM, M, WORK, LWORK,
     $                   IWORK, INFO )
*
*  -- LAPACK routine (version 3.3.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- April 2011                                                      --
*
*     .. Scalar Arguments ..
      CHARACTER          HOWMNY, JOB
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      INTEGER            IWORK( * )
      REAL               DIF( * ), S( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), VL( LDVL, * ),
     $                   VR( LDVR, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CTGSNA estimates reciprocal condition numbers for specified
*  eigenvalues and/or eigenvectors of a matrix pair (A, B).
*
*  (A, B) must be in generalized Schur canonical form, that is, A and
*  B are both upper triangular.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          Specifies whether condition numbers are required for
*          eigenvalues (S) or eigenvectors (DIF):
*          = 'E': for eigenvalues only (S);
*          = 'V': for eigenvectors only (DIF);
*          = 'B': for both eigenvalues and eigenvectors (S and DIF).
*
*  HOWMNY  (input) CHARACTER*1
*          = 'A': compute condition numbers for all eigenpairs;
*          = 'S': compute condition numbers for selected eigenpairs
*                 specified by the array SELECT.
*
*  SELECT  (input) LOGICAL array, dimension (N)
*          If HOWMNY = 'S', SELECT specifies the eigenpairs for which
*          condition numbers are required. To select condition numbers
*          for the corresponding j-th eigenvalue and/or eigenvector,
*          SELECT(j) must be set to .TRUE..
*          If HOWMNY = 'A', SELECT is not referenced.
*
*  N       (input) INTEGER
*          The order of the square matrix pair (A, B). N >= 0.
*
*  A       (input) COMPLEX array, dimension (LDA,N)
*          The upper triangular matrix A in the pair (A,B).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,N).
*
*  B       (input) COMPLEX array, dimension (LDB,N)
*          The upper triangular matrix B in the pair (A, B).
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,N).
*
*  VL      (input) COMPLEX array, dimension (LDVL,M)
*          IF JOB = 'E' or 'B', VL must contain left eigenvectors of
*          (A, B), corresponding to the eigenpairs specified by HOWMNY
*          and SELECT.  The eigenvectors must be stored in consecutive
*          columns of VL, as returned by CTGEVC.
*          If JOB = 'V', VL is not referenced.
*
*  LDVL    (input) INTEGER
*          The leading dimension of the array VL. LDVL >= 1; and
*          If JOB = 'E' or 'B', LDVL >= N.
*
*  VR      (input) COMPLEX array, dimension (LDVR,M)
*          IF JOB = 'E' or 'B', VR must contain right eigenvectors of
*          (A, B), corresponding to the eigenpairs specified by HOWMNY
*          and SELECT.  The eigenvectors must be stored in consecutive
*          columns of VR, as returned by CTGEVC.
*          If JOB = 'V', VR is not referenced.
*
*  LDVR    (input) INTEGER
*          The leading dimension of the array VR. LDVR >= 1;
*          If JOB = 'E' or 'B', LDVR >= N.
*
*  S       (output) REAL array, dimension (MM)
*          If JOB = 'E' or 'B', the reciprocal condition numbers of the
*          selected eigenvalues, stored in consecutive elements of the
*          array.
*          If JOB = 'V', S is not referenced.
*
*  DIF     (output) REAL array, dimension (MM)
*          If JOB = 'V' or 'B', the estimated reciprocal condition
*          numbers of the selected eigenvectors, stored in consecutive
*          elements of the array.
*          If the eigenvalues cannot be reordered to compute DIF(j),
*          DIF(j) is set to 0; this can only occur when the true value
*          would be very small anyway.
*          For each eigenvalue/vector specified by SELECT, DIF stores
*          a Frobenius norm-based estimate of Difl.
*          If JOB = 'E', DIF is not referenced.
*
*  MM      (input) INTEGER
*          The number of elements in the arrays S and DIF. MM >= M.
*
*  M       (output) INTEGER
*          The number of elements of the arrays S and DIF used to store
*          the specified condition numbers; for each selected eigenvalue
*          one element is used. If HOWMNY = 'A', M is set to N.
*
*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK  (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,N).
*          If JOB = 'V' or 'B', LWORK >= max(1,2*N*N).
*
*  IWORK   (workspace) INTEGER array, dimension (N+2)
*          If JOB = 'E', IWORK is not referenced.
*
*  INFO    (output) INTEGER
*          = 0: Successful exit
*          < 0: If INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The reciprocal of the condition number of the i-th generalized
*  eigenvalue w = (a, b) is defined as
*
*          S(I) = (|v**HAu|**2 + |v**HBu|**2)**(1/2) / (norm(u)*norm(v))
*
*  where u and v are the right and left eigenvectors of (A, B)
*  corresponding to w; |z| denotes the absolute value of the complex
*  number, and norm(u) denotes the 2-norm of the vector u. The pair
*  (a, b) corresponds to an eigenvalue w = a/b (= v**HAu/v**HBu) of the
*  matrix pair (A, B). If both a and b equal zero, then (A,B) is
*  singular and S(I) = -1 is returned.
*
*  An approximate error bound on the chordal distance between the i-th
*  computed generalized eigenvalue w and the corresponding exact
*  eigenvalue lambda is
*
*          chord(w, lambda) <=   EPS * norm(A, B) / S(I),
*
*  where EPS is the machine precision.
*
*  The reciprocal of the condition number of the right eigenvector u
*  and left eigenvector v corresponding to the generalized eigenvalue w
*  is defined as follows. Suppose
*
*                   (A, B) = ( a   *  ) ( b  *  )  1
*                            ( 0  A22 ),( 0 B22 )  n-1
*                              1  n-1     1 n-1
*
*  Then the reciprocal condition number DIF(I) is
*
*          Difl[(a, b), (A22, B22)]  = sigma-min( Zl )
*
*  where sigma-min(Zl) denotes the smallest singular value of
*
*         Zl = [ kron(a, In-1) -kron(1, A22) ]
*              [ kron(b, In-1) -kron(1, B22) ].
*
*  Here In-1 is the identity matrix of size n-1 and X**H is the conjugate
*  transpose of X. kron(X, Y) is the Kronecker product between the
*  matrices X and Y.
*
*  We approximate the smallest singular value of Zl with an upper
*  bound. This is done by CLATDF.
*
*  An approximate error bound for a computed eigenvector VL(i) or
*  VR(i) is given by
*
*                      EPS * norm(A, B) / DIF(i).
*
*  See ref. [2-3] for more details and further references.
*
*  Based on contributions by
*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
*     Umea University, S-901 87 Umea, Sweden.
*
*  References
*  ==========
*
*  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the
*      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in
*      M.S. Moonen et al (eds), Linear Algebra for Large Scale and
*      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.
*
*  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified
*      Eigenvalues of a Regular Matrix Pair (A, B) and Condition
*      Estimation: Theory, Algorithms and Software, Report
*      UMINF - 94.04, Department of Computing Science, Umea University,
*      S-901 87 Umea, Sweden, 1994. Also as LAPACK Working Note 87.
*      To appear in Numerical Algorithms, 1996.
*
*  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software
*      for Solving the Generalized Sylvester Equation and Estimating the
*      Separation between Regular Matrix Pairs, Report UMINF - 93.23,
*      Department of Computing Science, Umea University, S-901 87 Umea,
*      Sweden, December 1993, Revised April 1994, Also as LAPACK Working
*      Note 75.
*      To appear in ACM Trans. on Math. Software, Vol 22, No 1, 1996.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      INTEGER            IDIFJB
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, IDIFJB = 3 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, SOMCON, WANTBH, WANTDF, WANTS
      INTEGER            I, IERR, IFST, ILST, K, KS, LWMIN, N1, N2
      REAL               BIGNUM, COND, EPS, LNRM, RNRM, SCALE, SMLNUM
      COMPLEX            YHAX, YHBX
*     ..
*     .. Local Arrays ..
      COMPLEX            DUMMY( 1 ), DUMMY1( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SCNRM2, SLAMCH, SLAPY2
      COMPLEX            CDOTC
      EXTERNAL           LSAME, SCNRM2, SLAMCH, SLAPY2, CDOTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMV, CLACPY, CTGEXC, CTGSYL, SLABAD, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CMPLX, MAX
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      WANTBH = LSAME( JOB, 'B' )
      WANTS = LSAME( JOB, 'E' ) .OR. WANTBH
      WANTDF = LSAME( JOB, 'V' ) .OR. WANTBH
*
      SOMCON = LSAME( HOWMNY, 'S' )
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
*
      IF( .NOT.WANTS .AND. .NOT.WANTDF ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( HOWMNY, 'A' ) .AND. .NOT.SOMCON ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( WANTS .AND. LDVL.LT.N ) THEN
         INFO = -10
      ELSE IF( WANTS .AND. LDVR.LT.N ) THEN
         INFO = -12
      ELSE
*
*        Set M to the number of eigenpairs for which condition numbers
*        are required, and test MM.
*
         IF( SOMCON ) THEN
            M = 0
            DO 10 K = 1, N
               IF( SELECT( K ) )
     $            M = M + 1
   10       CONTINUE
         ELSE
            M = N
         END IF
*
         IF( N.EQ.0 ) THEN
            LWMIN = 1
         ELSE IF( LSAME( JOB, 'V' ) .OR. LSAME( JOB, 'B' ) ) THEN
            LWMIN = 2*N*N
         ELSE
            LWMIN = N
         END IF
         WORK( 1 ) = LWMIN
*
         IF( MM.LT.M ) THEN
            INFO = -15
         ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -18
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTGSNA', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Get machine constants
*
      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM
      CALL SLABAD( SMLNUM, BIGNUM )
      KS = 0
      DO 20 K = 1, N
*
*        Determine whether condition numbers are required for the k-th
*        eigenpair.
*
         IF( SOMCON ) THEN
            IF( .NOT.SELECT( K ) )
     $         GO TO 20
         END IF
*
         KS = KS + 1
*
         IF( WANTS ) THEN
*
*           Compute the reciprocal condition number of the k-th
*           eigenvalue.
*
            RNRM = SCNRM2( N, VR( 1, KS ), 1 )
            LNRM = SCNRM2( N, VL( 1, KS ), 1 )
            CALL CGEMV( 'N', N, N, CMPLX( ONE, ZERO ), A, LDA,
     $                  VR( 1, KS ), 1, CMPLX( ZERO, ZERO ), WORK, 1 )
            YHAX = CDOTC( N, WORK, 1, VL( 1, KS ), 1 )
            CALL CGEMV( 'N', N, N, CMPLX( ONE, ZERO ), B, LDB,
     $                  VR( 1, KS ), 1, CMPLX( ZERO, ZERO ), WORK, 1 )
            YHBX = CDOTC( N, WORK, 1, VL( 1, KS ), 1 )
            COND = SLAPY2( ABS( YHAX ), ABS( YHBX ) )
            IF( COND.EQ.ZERO ) THEN
               S( KS ) = -ONE
            ELSE
               S( KS ) = COND / ( RNRM*LNRM )
            END IF
         END IF
*
         IF( WANTDF ) THEN
            IF( N.EQ.1 ) THEN
               DIF( KS ) = SLAPY2( ABS( A( 1, 1 ) ), ABS( B( 1, 1 ) ) )
            ELSE
*
*              Estimate the reciprocal condition number of the k-th
*              eigenvectors.
*
*              Copy the matrix (A, B) to the array WORK and move the
*              (k,k)th pair to the (1,1) position.
*
               CALL CLACPY( 'Full', N, N, A, LDA, WORK, N )
               CALL CLACPY( 'Full', N, N, B, LDB, WORK( N*N+1 ), N )
               IFST = K
               ILST = 1
*
               CALL CTGEXC( .FALSE., .FALSE., N, WORK, N, WORK( N*N+1 ),
     $                      N, DUMMY, 1, DUMMY1, 1, IFST, ILST, IERR )
*
               IF( IERR.GT.0 ) THEN
*
*                 Ill-conditioned problem - swap rejected.
*
                  DIF( KS ) = ZERO
               ELSE
*
*                 Reordering successful, solve generalized Sylvester
*                 equation for R and L,
*                            A22 * R - L * A11 = A12
*                            B22 * R - L * B11 = B12,
*                 and compute estimate of Difl[(A11,B11), (A22, B22)].
*
                  N1 = 1
                  N2 = N - N1
                  I = N*N + 1
                  CALL CTGSYL( 'N', IDIFJB, N2, N1, WORK( N*N1+N1+1 ),
     $                         N, WORK, N, WORK( N1+1 ), N,
     $                         WORK( N*N1+N1+I ), N, WORK( I ), N,
     $                         WORK( N1+I ), N, SCALE, DIF( KS ), DUMMY,
     $                         1, IWORK, IERR )
               END IF
            END IF
         END IF
*
   20 CONTINUE
      WORK( 1 ) = LWMIN
      RETURN
*
*     End of CTGSNA
*
      END
