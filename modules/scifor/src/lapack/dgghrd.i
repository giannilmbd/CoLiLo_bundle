# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/dgghrd.f"
      SUBROUTINE DGGHRD( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q,
     $                   LDQ, Z, LDZ, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ
      INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
     $                   Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  DGGHRD reduces a pair of real matrices (A,B) to generalized upper
*  Hessenberg form using orthogonal transformations, where A is a
*  general matrix and B is upper triangular.  The form of the
*  generalized eigenvalue problem is
*     A*x = lambda*B*x,
*  and B is typically made upper triangular by computing its QR
*  factorization and moving the orthogonal matrix Q to the left side
*  of the equation.
*
*  This subroutine simultaneously reduces A to a Hessenberg matrix H:
*     Q**T*A*Z = H
*  and transforms B to another upper triangular matrix T:
*     Q**T*B*Z = T
*  in order to reduce the problem to its standard form
*     H*y = lambda*T*y
*  where y = Z**T*x.
*
*  The orthogonal matrices Q and Z are determined as products of Givens
*  rotations.  They may either be formed explicitly, or they may be
*  postmultiplied into input matrices Q1 and Z1, so that
*
*       Q1 * A * Z1**T = (Q1*Q) * H * (Z1*Z)**T
*
*       Q1 * B * Z1**T = (Q1*Q) * T * (Z1*Z)**T
*
*  If Q1 is the orthogonal matrix from the QR factorization of B in the
*  original equation A*x = lambda*B*x, then DGGHRD reduces the original
*  problem to generalized Hessenberg form.
*
*  Arguments
*  =========
*
*  COMPQ   (input) CHARACTER*1
*          = 'N': do not compute Q;
*          = 'I': Q is initialized to the unit matrix, and the
*                 orthogonal matrix Q is returned;
*          = 'V': Q must contain an orthogonal matrix Q1 on entry,
*                 and the product Q1*Q is returned.
*
*  COMPZ   (input) CHARACTER*1
*          = 'N': do not compute Z;
*          = 'I': Z is initialized to the unit matrix, and the
*                 orthogonal matrix Z is returned;
*          = 'V': Z must contain an orthogonal matrix Z1 on entry,
*                 and the product Z1*Z is returned.
*
*  N       (input) INTEGER
*          The order of the matrices A and B.  N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          ILO and IHI mark the rows and columns of A which are to be
*          reduced.  It is assumed that A is already upper triangular
*          in rows and columns 1:ILO-1 and IHI+1:N.  ILO and IHI are
*          normally set by a previous call to SGGBAL; otherwise they
*          should be set to 1 and N respectively.
*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
*          On entry, the N-by-N general matrix to be reduced.
*          On exit, the upper triangle and the first subdiagonal of A
*          are overwritten with the upper Hessenberg matrix H, and the
*          rest is set to zero.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
*          On entry, the N-by-N upper triangular matrix B.
*          On exit, the upper triangular matrix T = Q**T B Z.  The
*          elements below the diagonal are set to zero.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ, N)
*          On entry, if COMPQ = 'V', the orthogonal matrix Q1,
*          typically from the QR factorization of B.
*          On exit, if COMPQ='I', the orthogonal matrix Q, and if
*          COMPQ = 'V', the product Q1*Q.
*          Not referenced if COMPQ='N'.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.
*          LDQ >= N if COMPQ='V' or 'I'; LDQ >= 1 otherwise.
*
*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)
*          On entry, if COMPZ = 'V', the orthogonal matrix Z1.
*          On exit, if COMPZ='I', the orthogonal matrix Z, and if
*          COMPZ = 'V', the product Z1*Z.
*          Not referenced if COMPZ='N'.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.
*          LDZ >= N if COMPZ='V' or 'I'; LDZ >= 1 otherwise.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  This routine reduces A to Hessenberg and B to triangular form by
*  an unblocked reduction, as described in _Matrix_Computations_,
*  by Golub and Van Loan (Johns Hopkins Press.)
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILQ, ILZ
      INTEGER            ICOMPQ, ICOMPZ, JCOL, JROW
      DOUBLE PRECISION   C, S, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARTG, DLASET, DROT, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Decode COMPQ
*
      IF( LSAME( COMPQ, 'N' ) ) THEN
         ILQ = .FALSE.
         ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'V' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 2
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 3
      ELSE
         ICOMPQ = 0
      END IF
*
*     Decode COMPZ
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ILZ = .FALSE.
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 2
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 3
      ELSE
         ICOMPZ = 0
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF( ICOMPQ.LE.0 ) THEN
         INFO = -1
      ELSE IF( ICOMPZ.LE.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -4
      ELSE IF( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( ( ILQ .AND. LDQ.LT.N ) .OR. LDQ.LT.1 ) THEN
         INFO = -11
      ELSE IF( ( ILZ .AND. LDZ.LT.N ) .OR. LDZ.LT.1 ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGHRD', -INFO )
         RETURN
      END IF
*
*     Initialize Q and Z if desired.
*
      IF( ICOMPQ.EQ.3 )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDQ )
      IF( ICOMPZ.EQ.3 )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
*     Zero out lower triangle of B
*
      DO 20 JCOL = 1, N - 1
         DO 10 JROW = JCOL + 1, N
            B( JROW, JCOL ) = ZERO
   10    CONTINUE
   20 CONTINUE
*
*     Reduce A and B
*
      DO 40 JCOL = ILO, IHI - 2
*
         DO 30 JROW = IHI, JCOL + 2, -1
*
*           Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL)
*
            TEMP = A( JROW-1, JCOL )
            CALL DLARTG( TEMP, A( JROW, JCOL ), C, S,
     $                   A( JROW-1, JCOL ) )
            A( JROW, JCOL ) = ZERO
            CALL DROT( N-JCOL, A( JROW-1, JCOL+1 ), LDA,
     $                 A( JROW, JCOL+1 ), LDA, C, S )
            CALL DROT( N+2-JROW, B( JROW-1, JROW-1 ), LDB,
     $                 B( JROW, JROW-1 ), LDB, C, S )
            IF( ILQ )
     $         CALL DROT( N, Q( 1, JROW-1 ), 1, Q( 1, JROW ), 1, C, S )
*
*           Step 2: rotate columns JROW, JROW-1 to kill B(JROW,JROW-1)
*
            TEMP = B( JROW, JROW )
            CALL DLARTG( TEMP, B( JROW, JROW-1 ), C, S,
     $                   B( JROW, JROW ) )
            B( JROW, JROW-1 ) = ZERO
            CALL DROT( IHI, A( 1, JROW ), 1, A( 1, JROW-1 ), 1, C, S )
            CALL DROT( JROW-1, B( 1, JROW ), 1, B( 1, JROW-1 ), 1, C,
     $                 S )
            IF( ILZ )
     $         CALL DROT( N, Z( 1, JROW ), 1, Z( 1, JROW-1 ), 1, C, S )
   30    CONTINUE
   40 CONTINUE
*
      RETURN
*
*     End of DGGHRD
*
      END
