# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/chbevx.f"
      SUBROUTINE CHBEVX( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, Q, LDQ, VL,
     $                   VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK,
     $                   IWORK, IFAIL, INFO )
*
*  -- LAPACK driver routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N
      REAL               ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            IFAIL( * ), IWORK( * )
      REAL               RWORK( * ), W( * )
      COMPLEX            AB( LDAB, * ), Q( LDQ, * ), WORK( * ),
     $                   Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  CHBEVX computes selected eigenvalues and, optionally, eigenvectors
*  of a complex Hermitian band matrix A.  Eigenvalues and eigenvectors
*  can be selected by specifying either a range of values or a range of
*  indices for the desired eigenvalues.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  RANGE   (input) CHARACTER*1
*          = 'A': all eigenvalues will be found;
*          = 'V': all eigenvalues in the half-open interval (VL,VU]
*                 will be found;
*          = 'I': the IL-th through IU-th eigenvalues will be found.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
*
*  AB      (input/output) COMPLEX array, dimension (LDAB, N)
*          On entry, the upper or lower triangle of the Hermitian band
*          matrix A, stored in the first KD+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*
*          On exit, AB is overwritten by values generated during the
*          reduction to tridiagonal form.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD + 1.
*
*  Q       (output) COMPLEX array, dimension (LDQ, N)
*          If JOBZ = 'V', the N-by-N unitary matrix used in the
*                          reduction to tridiagonal form.
*          If JOBZ = 'N', the array Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.  If JOBZ = 'V', then
*          LDQ >= max(1,N).
*
*  VL      (input) REAL
*  VU      (input) REAL
*          If RANGE='V', the lower and upper bounds of the interval to
*          be searched for eigenvalues. VL < VU.
*          Not referenced if RANGE = 'A' or 'I'.
*
*  IL      (input) INTEGER
*  IU      (input) INTEGER
*          If RANGE='I', the indices (in ascending order) of the
*          smallest and largest eigenvalues to be returned.
*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
*          Not referenced if RANGE = 'A' or 'V'.
*
*  ABSTOL  (input) REAL
*          The absolute error tolerance for the eigenvalues.
*          An approximate eigenvalue is accepted as converged
*          when it is determined to lie in an interval [a,b]
*          of width less than or equal to
*
*                  ABSTOL + EPS *   max( |a|,|b| ) ,
*
*          where EPS is the machine precision.  If ABSTOL is less than
*          or equal to zero, then  EPS*|T|  will be used in its place,
*          where |T| is the 1-norm of the tridiagonal matrix obtained
*          by reducing AB to tridiagonal form.
*
*          Eigenvalues will be computed most accurately when ABSTOL is
*          set to twice the underflow threshold 2*SLAMCH('S'), not zero.
*          If this routine returns with INFO>0, indicating that some
*          eigenvectors did not converge, try setting ABSTOL to
*          2*SLAMCH('S').
*
*          See "Computing Small Singular Values of Bidiagonal Matrices
*          with Guaranteed High Relative Accuracy," by Demmel and
*          Kahan, LAPACK Working Note #3.
*
*  M       (output) INTEGER
*          The total number of eigenvalues found.  0 <= M <= N.
*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
*
*  W       (output) REAL array, dimension (N)
*          The first M elements contain the selected eigenvalues in
*          ascending order.
*
*  Z       (output) COMPLEX array, dimension (LDZ, max(1,M))
*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
*          contain the orthonormal eigenvectors of the matrix A
*          corresponding to the selected eigenvalues, with the i-th
*          column of Z holding the eigenvector associated with W(i).
*          If an eigenvector fails to converge, then that column of Z
*          contains the latest approximation to the eigenvector, and the
*          index of the eigenvector is returned in IFAIL.
*          If JOBZ = 'N', then Z is not referenced.
*          Note: the user must ensure that at least max(1,M) columns are
*          supplied in the array Z; if RANGE = 'V', the exact value of M
*          is not known in advance and an upper bound must be used.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  WORK    (workspace) COMPLEX array, dimension (N)
*
*  RWORK   (workspace) REAL array, dimension (7*N)
*
*  IWORK   (workspace) INTEGER array, dimension (5*N)
*
*  IFAIL   (output) INTEGER array, dimension (N)
*          If JOBZ = 'V', then if INFO = 0, the first M elements of
*          IFAIL are zero.  If INFO > 0, then IFAIL contains the
*          indices of the eigenvectors that failed to converge.
*          If JOBZ = 'N', then IFAIL is not referenced.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, then i eigenvectors failed to converge.
*                Their indices are stored in array IFAIL.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E0, 0.0E0 ),
     $                   CONE = ( 1.0E0, 0.0E0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, LOWER, TEST, VALEIG, WANTZ
      CHARACTER          ORDER
      INTEGER            I, IINFO, IMAX, INDD, INDE, INDEE, INDIBL,
     $                   INDISP, INDIWK, INDRWK, INDWRK, ISCALE, ITMP1,
     $                   J, JJ, NSPLIT
      REAL               ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN,
     $                   SIGMA, SMLNUM, TMP1, VLL, VUU
      COMPLEX            CTMP1
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               CLANHB, SLAMCH
      EXTERNAL           LSAME, CLANHB, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CCOPY, CGEMV, CHBTRD, CLACPY, CLASCL, CSTEIN,
     $                   CSTEQR, CSWAP, SCOPY, SSCAL, SSTEBZ, SSTERF,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
      LOWER = LSAME( UPLO, 'L' )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( KD.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -7
      ELSE IF( WANTZ .AND. LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE
         IF( VALEIG ) THEN
            IF( N.GT.0 .AND. VU.LE.VL )
     $         INFO = -11
         ELSE IF( INDEIG ) THEN
            IF( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) THEN
               INFO = -12
            ELSE IF( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) THEN
               INFO = -13
            END IF
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) )
     $     INFO = -18
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHBEVX', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      M = 0
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         M = 1
         IF( LOWER ) THEN
            CTMP1 = AB( 1, 1 )
         ELSE
            CTMP1 = AB( KD+1, 1 )
         END IF
         TMP1 = REAL( CTMP1 )
         IF( VALEIG ) THEN
            IF( .NOT.( VL.LT.TMP1 .AND. VU.GE.TMP1 ) )
     $         M = 0
         END IF
         IF( M.EQ.1 ) THEN
            W( 1 ) = CTMP1
            IF( WANTZ )
     $         Z( 1, 1 ) = CONE
         END IF
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = SLAMCH( 'Safe minimum' )
      EPS = SLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
*
*     Scale matrix to allowable range, if necessary.
*
      ISCALE = 0
      ABSTLL = ABSTOL
      IF ( VALEIG ) THEN
         VLL = VL
         VUU = VU
      ELSE
         VLL = ZERO
         VUU = ZERO
      ENDIF
      ANRM = CLANHB( 'M', UPLO, N, KD, AB, LDAB, RWORK )
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         IF( LOWER ) THEN
            CALL CLASCL( 'B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO )
         ELSE
            CALL CLASCL( 'Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO )
         END IF
         IF( ABSTOL.GT.0 )
     $      ABSTLL = ABSTOL*SIGMA
         IF( VALEIG ) THEN
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         END IF
      END IF
*
*     Call CHBTRD to reduce Hermitian band matrix to tridiagonal form.
*
      INDD = 1
      INDE = INDD + N
      INDRWK = INDE + N
      INDWRK = 1
      CALL CHBTRD( JOBZ, UPLO, N, KD, AB, LDAB, RWORK( INDD ),
     $             RWORK( INDE ), Q, LDQ, WORK( INDWRK ), IINFO )
*
*     If all eigenvalues are desired and ABSTOL is less than or equal
*     to zero, then call SSTERF or CSTEQR.  If this fails for some
*     eigenvalue, then try SSTEBZ.
*
      TEST = .FALSE.
      IF (INDEIG) THEN
         IF (IL.EQ.1 .AND. IU.EQ.N) THEN
            TEST = .TRUE.
         END IF
      END IF
      IF ((ALLEIG .OR. TEST) .AND. (ABSTOL.LE.ZERO)) THEN
         CALL SCOPY( N, RWORK( INDD ), 1, W, 1 )
         INDEE = INDRWK + 2*N
         IF( .NOT.WANTZ ) THEN
            CALL SCOPY( N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 )
            CALL SSTERF( N, W, RWORK( INDEE ), INFO )
         ELSE
            CALL CLACPY( 'A', N, N, Q, LDQ, Z, LDZ )
            CALL SCOPY( N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 )
            CALL CSTEQR( JOBZ, N, W, RWORK( INDEE ), Z, LDZ,
     $                   RWORK( INDRWK ), INFO )
            IF( INFO.EQ.0 ) THEN
               DO 10 I = 1, N
                  IFAIL( I ) = 0
   10          CONTINUE
            END IF
         END IF
         IF( INFO.EQ.0 ) THEN
            M = N
            GO TO 30
         END IF
         INFO = 0
      END IF
*
*     Otherwise, call SSTEBZ and, if eigenvectors are desired, CSTEIN.
*
      IF( WANTZ ) THEN
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      END IF
      INDIBL = 1
      INDISP = INDIBL + N
      INDIWK = INDISP + N
      CALL SSTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL,
     $             RWORK( INDD ), RWORK( INDE ), M, NSPLIT, W,
     $             IWORK( INDIBL ), IWORK( INDISP ), RWORK( INDRWK ),
     $             IWORK( INDIWK ), INFO )
*
      IF( WANTZ ) THEN
         CALL CSTEIN( N, RWORK( INDD ), RWORK( INDE ), M, W,
     $                IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ,
     $                RWORK( INDRWK ), IWORK( INDIWK ), IFAIL, INFO )
*
*        Apply unitary matrix used in reduction to tridiagonal
*        form to eigenvectors returned by CSTEIN.
*
         DO 20 J = 1, M
            CALL CCOPY( N, Z( 1, J ), 1, WORK( 1 ), 1 )
            CALL CGEMV( 'N', N, N, CONE, Q, LDQ, WORK, 1, CZERO,
     $                  Z( 1, J ), 1 )
   20    CONTINUE
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
   30 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = M
         ELSE
            IMAX = INFO - 1
         END IF
         CALL SSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
*
*     If eigenvalues are not in order, then sort them, along with
*     eigenvectors.
*
      IF( WANTZ ) THEN
         DO 50 J = 1, M - 1
            I = 0
            TMP1 = W( J )
            DO 40 JJ = J + 1, M
               IF( W( JJ ).LT.TMP1 ) THEN
                  I = JJ
                  TMP1 = W( JJ )
               END IF
   40       CONTINUE
*
            IF( I.NE.0 ) THEN
               ITMP1 = IWORK( INDIBL+I-1 )
               W( I ) = W( J )
               IWORK( INDIBL+I-1 ) = IWORK( INDIBL+J-1 )
               W( J ) = TMP1
               IWORK( INDIBL+J-1 ) = ITMP1
               CALL CSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
               IF( INFO.NE.0 ) THEN
                  ITMP1 = IFAIL( I )
                  IFAIL( I ) = IFAIL( J )
                  IFAIL( J ) = ITMP1
               END IF
            END IF
   50    CONTINUE
      END IF
*
      RETURN
*
*     End of CHBEVX
*
      END
