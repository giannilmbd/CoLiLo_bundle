# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/cggesx.f"
      SUBROUTINE CGGESX( JOBVSL, JOBVSR, SORT, SELCTG, SENSE, N, A, LDA,
     $                   B, LDB, SDIM, ALPHA, BETA, VSL, LDVSL, VSR,
     $                   LDVSR, RCONDE, RCONDV, WORK, LWORK, RWORK,
     $                   IWORK, LIWORK, BWORK, INFO )
*
*  -- LAPACK driver routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          JOBVSL, JOBVSR, SENSE, SORT
      INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LIWORK, LWORK, N,
     $                   SDIM
*     ..
*     .. Array Arguments ..
      LOGICAL            BWORK( * )
      INTEGER            IWORK( * )
      REAL               RCONDE( 2 ), RCONDV( 2 ), RWORK( * )
      COMPLEX            A( LDA, * ), ALPHA( * ), B( LDB, * ),
     $                   BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ),
     $                   WORK( * )
*     ..
*     .. Function Arguments ..
      LOGICAL            SELCTG
      EXTERNAL           SELCTG
*     ..
*
*  Purpose
*  =======
*
*  CGGESX computes for a pair of N-by-N complex nonsymmetric matrices
*  (A,B), the generalized eigenvalues, the complex Schur form (S,T),
*  and, optionally, the left and/or right matrices of Schur vectors (VSL
*  and VSR).  This gives the generalized Schur factorization
*
*       (A,B) = ( (VSL) S (VSR)**H, (VSL) T (VSR)**H )
*
*  where (VSR)**H is the conjugate-transpose of VSR.
*
*  Optionally, it also orders the eigenvalues so that a selected cluster
*  of eigenvalues appears in the leading diagonal blocks of the upper
*  triangular matrix S and the upper triangular matrix T; computes
*  a reciprocal condition number for the average of the selected
*  eigenvalues (RCONDE); and computes a reciprocal condition number for
*  the right and left deflating subspaces corresponding to the selected
*  eigenvalues (RCONDV). The leading columns of VSL and VSR then form
*  an orthonormal basis for the corresponding left and right eigenspaces
*  (deflating subspaces).
*
*  A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
*  or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
*  usually represented as the pair (alpha,beta), as there is a
*  reasonable interpretation for beta=0 or for both being zero.
*
*  A pair of matrices (S,T) is in generalized complex Schur form if T is
*  upper triangular with non-negative diagonal and S is upper
*  triangular.
*
*  Arguments
*  =========
*
*  JOBVSL  (input) CHARACTER*1
*          = 'N':  do not compute the left Schur vectors;
*          = 'V':  compute the left Schur vectors.
*
*  JOBVSR  (input) CHARACTER*1
*          = 'N':  do not compute the right Schur vectors;
*          = 'V':  compute the right Schur vectors.
*
*  SORT    (input) CHARACTER*1
*          Specifies whether or not to order the eigenvalues on the
*          diagonal of the generalized Schur form.
*          = 'N':  Eigenvalues are not ordered;
*          = 'S':  Eigenvalues are ordered (see SELCTG).
*
*  SELCTG  (external procedure) LOGICAL FUNCTION of two COMPLEX arguments
*          SELCTG must be declared EXTERNAL in the calling subroutine.
*          If SORT = 'N', SELCTG is not referenced.
*          If SORT = 'S', SELCTG is used to select eigenvalues to sort
*          to the top left of the Schur form.
*          Note that a selected complex eigenvalue may no longer satisfy
*          SELCTG(ALPHA(j),BETA(j)) = .TRUE. after ordering, since
*          ordering may change the value of complex eigenvalues
*          (especially if the eigenvalue is ill-conditioned), in this
*          case INFO is set to N+3 see INFO below).
*
*  SENSE   (input) CHARACTER*1
*          Determines which reciprocal condition numbers are computed.
*          = 'N' : None are computed;
*          = 'E' : Computed for average of selected eigenvalues only;
*          = 'V' : Computed for selected deflating subspaces only;
*          = 'B' : Computed for both.
*          If SENSE = 'E', 'V', or 'B', SORT must equal 'S'.
*
*  N       (input) INTEGER
*          The order of the matrices A, B, VSL, and VSR.  N >= 0.
*
*  A       (input/output) COMPLEX array, dimension (LDA, N)
*          On entry, the first of the pair of matrices.
*          On exit, A has been overwritten by its generalized Schur
*          form S.
*
*  LDA     (input) INTEGER
*          The leading dimension of A.  LDA >= max(1,N).
*
*  B       (input/output) COMPLEX array, dimension (LDB, N)
*          On entry, the second of the pair of matrices.
*          On exit, B has been overwritten by its generalized Schur
*          form T.
*
*  LDB     (input) INTEGER
*          The leading dimension of B.  LDB >= max(1,N).
*
*  SDIM    (output) INTEGER
*          If SORT = 'N', SDIM = 0.
*          If SORT = 'S', SDIM = number of eigenvalues (after sorting)
*          for which SELCTG is true.
*
*  ALPHA   (output) COMPLEX array, dimension (N)
*  BETA    (output) COMPLEX array, dimension (N)
*          On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the
*          generalized eigenvalues.  ALPHA(j) and BETA(j),j=1,...,N  are
*          the diagonals of the complex Schur form (S,T).  BETA(j) will
*          be non-negative real.
*
*          Note: the quotients ALPHA(j)/BETA(j) may easily over- or
*          underflow, and BETA(j) may even be zero.  Thus, the user
*          should avoid naively computing the ratio alpha/beta.
*          However, ALPHA will be always less than and usually
*          comparable with norm(A) in magnitude, and BETA always less
*          than and usually comparable with norm(B).
*
*  VSL     (output) COMPLEX array, dimension (LDVSL,N)
*          If JOBVSL = 'V', VSL will contain the left Schur vectors.
*          Not referenced if JOBVSL = 'N'.
*
*  LDVSL   (input) INTEGER
*          The leading dimension of the matrix VSL. LDVSL >=1, and
*          if JOBVSL = 'V', LDVSL >= N.
*
*  VSR     (output) COMPLEX array, dimension (LDVSR,N)
*          If JOBVSR = 'V', VSR will contain the right Schur vectors.
*          Not referenced if JOBVSR = 'N'.
*
*  LDVSR   (input) INTEGER
*          The leading dimension of the matrix VSR. LDVSR >= 1, and
*          if JOBVSR = 'V', LDVSR >= N.
*
*  RCONDE  (output) REAL array, dimension ( 2 )
*          If SENSE = 'E' or 'B', RCONDE(1) and RCONDE(2) contain the
*          reciprocal condition numbers for the average of the selected
*          eigenvalues.
*          Not referenced if SENSE = 'N' or 'V'.
*
*  RCONDV  (output) REAL array, dimension ( 2 )
*          If SENSE = 'V' or 'B', RCONDV(1) and RCONDV(2) contain the
*          reciprocal condition number for the selected deflating
*          subspaces.
*          Not referenced if SENSE = 'N' or 'E'.
*
*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If N = 0, LWORK >= 1, else if SENSE = 'E', 'V', or 'B',
*          LWORK >= MAX(1,2*N,2*SDIM*(N-SDIM)), else
*          LWORK >= MAX(1,2*N).  Note that 2*SDIM*(N-SDIM) <= N*N/2.
*          Note also that an error is only returned if
*          LWORK < MAX(1,2*N), but if SENSE = 'E' or 'V' or 'B' this may
*          not be large enough.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the bound on the optimal size of the WORK
*          array and the minimum size of the IWORK array, returns these
*          values as the first entries of the WORK and IWORK arrays, and
*          no error message related to LWORK or LIWORK is issued by
*          XERBLA.
*
*  RWORK   (workspace) REAL array, dimension ( 8*N )
*          Real workspace.
*
*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
*          On exit, if INFO = 0, IWORK(1) returns the minimum LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array WORK.
*          If SENSE = 'N' or N = 0, LIWORK >= 1, otherwise
*          LIWORK >= N+2.
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the bound on the optimal size of the
*          WORK array and the minimum size of the IWORK array, returns
*          these values as the first entries of the WORK and IWORK
*          arrays, and no error message related to LWORK or LIWORK is
*          issued by XERBLA.
*
*  BWORK   (workspace) LOGICAL array, dimension (N)
*          Not referenced if SORT = 'N'.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          = 1,...,N:
*                The QZ iteration failed.  (A,B) are not in Schur
*                form, but ALPHA(j) and BETA(j) should be correct for
*                j=INFO+1,...,N.
*          > N:  =N+1: other than QZ iteration failed in CHGEQZ
*                =N+2: after reordering, roundoff changed values of
*                      some complex eigenvalues so that leading
*                      eigenvalues in the Generalized Schur form no
*                      longer satisfy SELCTG=.TRUE.  This could also
*                      be caused due to scaling.
*                =N+3: reordering failed in CTGSEN.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ),
     $                   CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            CURSL, ILASCL, ILBSCL, ILVSL, ILVSR, LASTSL,
     $                   LQUERY, WANTSB, WANTSE, WANTSN, WANTST, WANTSV
      INTEGER            I, ICOLS, IERR, IHI, IJOB, IJOBVL, IJOBVR,
     $                   ILEFT, ILO, IRIGHT, IROWS, IRWRK, ITAU, IWRK,
     $                   LIWMIN, LWRK, MAXWRK, MINWRK
      REAL               ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, PL,
     $                   PR, SMLNUM
*     ..
*     .. Local Arrays ..
      REAL               DIF( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEQRF, CGGBAK, CGGBAL, CGGHRD, CHGEQZ, CLACPY,
     $                   CLASCL, CLASET, CTGSEN, CUNGQR, CUNMQR, SLABAD,
     $                   XERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      REAL               CLANGE, SLAMCH
      EXTERNAL           LSAME, ILAENV, CLANGE, SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Decode the input arguments
*
      IF( LSAME( JOBVSL, 'N' ) ) THEN
         IJOBVL = 1
         ILVSL = .FALSE.
      ELSE IF( LSAME( JOBVSL, 'V' ) ) THEN
         IJOBVL = 2
         ILVSL = .TRUE.
      ELSE
         IJOBVL = -1
         ILVSL = .FALSE.
      END IF
*
      IF( LSAME( JOBVSR, 'N' ) ) THEN
         IJOBVR = 1
         ILVSR = .FALSE.
      ELSE IF( LSAME( JOBVSR, 'V' ) ) THEN
         IJOBVR = 2
         ILVSR = .TRUE.
      ELSE
         IJOBVR = -1
         ILVSR = .FALSE.
      END IF
*
      WANTST = LSAME( SORT, 'S' )
      WANTSN = LSAME( SENSE, 'N' )
      WANTSE = LSAME( SENSE, 'E' )
      WANTSV = LSAME( SENSE, 'V' )
      WANTSB = LSAME( SENSE, 'B' )
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
      IF( WANTSN ) THEN
         IJOB = 0
      ELSE IF( WANTSE ) THEN
         IJOB = 1
      ELSE IF( WANTSV ) THEN
         IJOB = 2
      ELSE IF( WANTSB ) THEN
         IJOB = 4
      END IF
*
*     Test the input arguments
*
      INFO = 0
      IF( IJOBVL.LE.0 ) THEN
         INFO = -1
      ELSE IF( IJOBVR.LE.0 ) THEN
         INFO = -2
      ELSE IF( ( .NOT.WANTST ) .AND. ( .NOT.LSAME( SORT, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( .NOT.( WANTSN .OR. WANTSE .OR. WANTSV .OR. WANTSB ) .OR.
     $         ( .NOT.WANTST .AND. .NOT.WANTSN ) ) THEN
         INFO = -5
      ELSE IF( N.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDVSL.LT.1 .OR. ( ILVSL .AND. LDVSL.LT.N ) ) THEN
         INFO = -15
      ELSE IF( LDVSR.LT.1 .OR. ( ILVSR .AND. LDVSR.LT.N ) ) THEN
         INFO = -17
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       NB refers to the optimal block size for the immediately
*       following subroutine, as returned by ILAENV.)
*
      IF( INFO.EQ.0 ) THEN
         IF( N.GT.0) THEN
            MINWRK = 2*N
            MAXWRK = N*(1 + ILAENV( 1, 'CGEQRF', ' ', N, 1, N, 0 ) )
            MAXWRK = MAX( MAXWRK, N*( 1 +
     $                    ILAENV( 1, 'CUNMQR', ' ', N, 1, N, -1 ) ) )
            IF( ILVSL ) THEN
               MAXWRK = MAX( MAXWRK, N*( 1 +
     $                       ILAENV( 1, 'CUNGQR', ' ', N, 1, N, -1 ) ) )
            END IF
            LWRK = MAXWRK
            IF( IJOB.GE.1 )
     $         LWRK = MAX( LWRK, N*N/2 )
         ELSE
            MINWRK = 1
            MAXWRK = 1
            LWRK   = 1
         END IF
         WORK( 1 ) = LWRK
         IF( WANTSN .OR. N.EQ.0 ) THEN
            LIWMIN = 1
         ELSE
            LIWMIN = N + 2
         END IF
         IWORK( 1 ) = LIWMIN
*
         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -21
         ELSE IF( LIWORK.LT.LIWMIN  .AND. .NOT.LQUERY) THEN
            INFO = -24
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGGESX', -INFO )
         RETURN
      ELSE IF (LQUERY) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         SDIM = 0
         RETURN
      END IF
*
*     Get machine constants
*
      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL SLABAD( SMLNUM, BIGNUM )
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      ANRM = CLANGE( 'M', N, N, A, LDA, RWORK )
      ILASCL = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         ANRMTO = SMLNUM
         ILASCL = .TRUE.
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         ANRMTO = BIGNUM
         ILASCL = .TRUE.
      END IF
      IF( ILASCL )
     $   CALL CLASCL( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR )
*
*     Scale B if max element outside range [SMLNUM,BIGNUM]
*
      BNRM = CLANGE( 'M', N, N, B, LDB, RWORK )
      ILBSCL = .FALSE.
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
         BNRMTO = SMLNUM
         ILBSCL = .TRUE.
      ELSE IF( BNRM.GT.BIGNUM ) THEN
         BNRMTO = BIGNUM
         ILBSCL = .TRUE.
      END IF
      IF( ILBSCL )
     $   CALL CLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR )
*
*     Permute the matrix to make it more nearly triangular
*     (Real Workspace: need 6*N)
*
      ILEFT = 1
      IRIGHT = N + 1
      IRWRK = IRIGHT + N
      CALL CGGBAL( 'P', N, A, LDA, B, LDB, ILO, IHI, RWORK( ILEFT ),
     $             RWORK( IRIGHT ), RWORK( IRWRK ), IERR )
*
*     Reduce B to triangular form (QR decomposition of B)
*     (Complex Workspace: need N, prefer N*NB)
*
      IROWS = IHI + 1 - ILO
      ICOLS = N + 1 - ILO
      ITAU = 1
      IWRK = ITAU + IROWS
      CALL CGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ),
     $             WORK( IWRK ), LWORK+1-IWRK, IERR )
*
*     Apply the unitary transformation to matrix A
*     (Complex Workspace: need N, prefer N*NB)
*
      CALL CUNMQR( 'L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB,
     $             WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ),
     $             LWORK+1-IWRK, IERR )
*
*     Initialize VSL
*     (Complex Workspace: need N, prefer N*NB)
*
      IF( ILVSL ) THEN
         CALL CLASET( 'Full', N, N, CZERO, CONE, VSL, LDVSL )
         IF( IROWS.GT.1 ) THEN
            CALL CLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB,
     $                   VSL( ILO+1, ILO ), LDVSL )
         END IF
         CALL CUNGQR( IROWS, IROWS, IROWS, VSL( ILO, ILO ), LDVSL,
     $                WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )
      END IF
*
*     Initialize VSR
*
      IF( ILVSR )
     $   CALL CLASET( 'Full', N, N, CZERO, CONE, VSR, LDVSR )
*
*     Reduce to generalized Hessenberg form
*     (Workspace: none needed)
*
      CALL CGGHRD( JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, VSL,
     $             LDVSL, VSR, LDVSR, IERR )
*
      SDIM = 0
*
*     Perform QZ algorithm, computing Schur vectors if desired
*     (Complex Workspace: need N)
*     (Real Workspace:    need N)
*
      IWRK = ITAU
      CALL CHGEQZ( 'S', JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB,
     $             ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK( IWRK ),
     $             LWORK+1-IWRK, RWORK( IRWRK ), IERR )
      IF( IERR.NE.0 ) THEN
         IF( IERR.GT.0 .AND. IERR.LE.N ) THEN
            INFO = IERR
         ELSE IF( IERR.GT.N .AND. IERR.LE.2*N ) THEN
            INFO = IERR - N
         ELSE
            INFO = N + 1
         END IF
         GO TO 40
      END IF
*
*     Sort eigenvalues ALPHA/BETA and compute the reciprocal of
*     condition number(s)
*
      IF( WANTST ) THEN
*
*        Undo scaling on eigenvalues before SELCTGing
*
         IF( ILASCL )
     $      CALL CLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA, N, IERR )
         IF( ILBSCL )
     $      CALL CLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )
*
*        Select eigenvalues
*
         DO 10 I = 1, N
            BWORK( I ) = SELCTG( ALPHA( I ), BETA( I ) )
   10    CONTINUE
*
*        Reorder eigenvalues, transform Generalized Schur vectors, and
*        compute reciprocal condition numbers
*        (Complex Workspace: If IJOB >= 1, need MAX(1, 2*SDIM*(N-SDIM))
*                            otherwise, need 1 )
*
         CALL CTGSEN( IJOB, ILVSL, ILVSR, BWORK, N, A, LDA, B, LDB,
     $                ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, SDIM, PL, PR,
     $                DIF, WORK( IWRK ), LWORK-IWRK+1, IWORK, LIWORK,
     $                IERR )
*
         IF( IJOB.GE.1 )
     $      MAXWRK = MAX( MAXWRK, 2*SDIM*( N-SDIM ) )
         IF( IERR.EQ.-21 ) THEN
*
*            not enough complex workspace
*
            INFO = -21
         ELSE
            IF( IJOB.EQ.1 .OR. IJOB.EQ.4 ) THEN
               RCONDE( 1 ) = PL
               RCONDE( 2 ) = PR
            END IF
            IF( IJOB.EQ.2 .OR. IJOB.EQ.4 ) THEN
               RCONDV( 1 ) = DIF( 1 )
               RCONDV( 2 ) = DIF( 2 )
            END IF
            IF( IERR.EQ.1 )
     $         INFO = N + 3
         END IF
*
      END IF
*
*     Apply permutation to VSL and VSR
*     (Workspace: none needed)
*
      IF( ILVSL )
     $   CALL CGGBAK( 'P', 'L', N, ILO, IHI, RWORK( ILEFT ),
     $                RWORK( IRIGHT ), N, VSL, LDVSL, IERR )
*
      IF( ILVSR )
     $   CALL CGGBAK( 'P', 'R', N, ILO, IHI, RWORK( ILEFT ),
     $                RWORK( IRIGHT ), N, VSR, LDVSR, IERR )
*
*     Undo scaling
*
      IF( ILASCL ) THEN
         CALL CLASCL( 'U', 0, 0, ANRMTO, ANRM, N, N, A, LDA, IERR )
         CALL CLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA, N, IERR )
      END IF
*
      IF( ILBSCL ) THEN
         CALL CLASCL( 'U', 0, 0, BNRMTO, BNRM, N, N, B, LDB, IERR )
         CALL CLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )
      END IF
*
      IF( WANTST ) THEN
*
*        Check if reordering is correct
*
         LASTSL = .TRUE.
         SDIM = 0
         DO 30 I = 1, N
            CURSL = SELCTG( ALPHA( I ), BETA( I ) )
            IF( CURSL )
     $         SDIM = SDIM + 1
            IF( CURSL .AND. .NOT.LASTSL )
     $         INFO = N + 2
            LASTSL = CURSL
   30    CONTINUE
*
      END IF
*
   40 CONTINUE
*
      WORK( 1 ) = MAXWRK
      IWORK( 1 ) = LIWMIN
*
      RETURN
*
*     End of CGGESX
*
      END
