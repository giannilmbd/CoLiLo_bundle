# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/dspev.f"
      SUBROUTINE DSPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO )
*
*  -- LAPACK driver routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDZ, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   AP( * ), W( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  DSPEV computes all the eigenvalues and, optionally, eigenvectors of a
*  real symmetric matrix A in packed storage.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
*          On entry, the upper or lower triangle of the symmetric matrix
*          A, packed columnwise in a linear array.  The j-th column of A
*          is stored in the array AP as follows:
*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
*
*          On exit, AP is overwritten by values generated during the
*          reduction to tridiagonal form.  If UPLO = 'U', the diagonal
*          and first superdiagonal of the tridiagonal matrix T overwrite
*          the corresponding elements of A, and if UPLO = 'L', the
*          diagonal and first subdiagonal of T overwrite the
*          corresponding elements of A.
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)
*          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
*          eigenvectors of the matrix A, with the i-th column of Z
*          holding the eigenvector associated with W(i).
*          If JOBZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of an intermediate tridiagonal
*                form did not converge to zero.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            WANTZ
      INTEGER            IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANSP
      EXTERNAL           LSAME, DLAMCH, DLANSP
*     ..
*     .. External Subroutines ..
      EXTERNAL           DOPGTR, DSCAL, DSPTRD, DSTEQR, DSTERF, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LSAME( UPLO, 'U' ) .OR. LSAME( UPLO, 'L' ) ) )
     $          THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -7
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSPEV ', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         W( 1 ) = AP( 1 )
         IF( WANTZ )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
*
*     Scale matrix to allowable range, if necessary.
*
      ANRM = DLANSP( 'M', UPLO, N, AP, WORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         CALL DSCAL( ( N*( N+1 ) ) / 2, SIGMA, AP, 1 )
      END IF
*
*     Call DSPTRD to reduce symmetric packed matrix to tridiagonal form.
*
      INDE = 1
      INDTAU = INDE + N
      CALL DSPTRD( UPLO, N, AP, W, WORK( INDE ), WORK( INDTAU ), IINFO )
*
*     For eigenvalues only, call DSTERF.  For eigenvectors, first call
*     DOPGTR to generate the orthogonal matrix, then call DSTEQR.
*
      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, W, WORK( INDE ), INFO )
      ELSE
         INDWRK = INDTAU + N
         CALL DOPGTR( UPLO, N, AP, WORK( INDTAU ), Z, LDZ,
     $                WORK( INDWRK ), IINFO )
         CALL DSTEQR( JOBZ, N, W, WORK( INDE ), Z, LDZ, WORK( INDTAU ),
     $                INFO )
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
*
      RETURN
*
*     End of DSPEV
*
      END
