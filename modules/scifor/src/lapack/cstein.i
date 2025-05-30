# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/cstein.f"
      SUBROUTINE CSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK,
     $                   IWORK, IFAIL, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDZ, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IBLOCK( * ), IFAIL( * ), ISPLIT( * ),
     $                   IWORK( * )
      REAL               D( * ), E( * ), W( * ), WORK( * )
      COMPLEX            Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  CSTEIN computes the eigenvectors of a real symmetric tridiagonal
*  matrix T corresponding to specified eigenvalues, using inverse
*  iteration.
*
*  The maximum number of iterations allowed for each eigenvector is
*  specified by an internal parameter MAXITS (currently set to 5).
*
*  Although the eigenvectors are real, they are stored in a complex
*  array, which may be passed to CUNMTR or CUPMTR for back
*  transformation to the eigenvectors of a complex Hermitian matrix
*  which was reduced to tridiagonal form.
*
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input) REAL array, dimension (N)
*          The n diagonal elements of the tridiagonal matrix T.
*
*  E       (input) REAL array, dimension (N-1)
*          The (n-1) subdiagonal elements of the tridiagonal matrix
*          T, stored in elements 1 to N-1.
*
*  M       (input) INTEGER
*          The number of eigenvectors to be found.  0 <= M <= N.
*
*  W       (input) REAL array, dimension (N)
*          The first M elements of W contain the eigenvalues for
*          which eigenvectors are to be computed.  The eigenvalues
*          should be grouped by split-off block and ordered from
*          smallest to largest within the block.  ( The output array
*          W from SSTEBZ with ORDER = 'B' is expected here. )
*
*  IBLOCK  (input) INTEGER array, dimension (N)
*          The submatrix indices associated with the corresponding
*          eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to
*          the first submatrix from the top, =2 if W(i) belongs to
*          the second submatrix, etc.  ( The output array IBLOCK
*          from SSTEBZ is expected here. )
*
*  ISPLIT  (input) INTEGER array, dimension (N)
*          The splitting points, at which T breaks up into submatrices.
*          The first submatrix consists of rows/columns 1 to
*          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1
*          through ISPLIT( 2 ), etc.
*          ( The output array ISPLIT from SSTEBZ is expected here. )
*
*  Z       (output) COMPLEX array, dimension (LDZ, M)
*          The computed eigenvectors.  The eigenvector associated
*          with the eigenvalue W(i) is stored in the i-th column of
*          Z.  Any vector which fails to converge is set to its current
*          iterate after MAXITS iterations.
*          The imaginary parts of the eigenvectors are set to zero.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= max(1,N).
*
*  WORK    (workspace) REAL array, dimension (5*N)
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*
*  IFAIL   (output) INTEGER array, dimension (M)
*          On normal exit, all elements of IFAIL are zero.
*          If one or more eigenvectors fail to converge after
*          MAXITS iterations, then their indices are stored in
*          array IFAIL.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = i, then i eigenvectors failed to converge
*               in MAXITS iterations.  Their indices are stored in
*               array IFAIL.
*
*  Internal Parameters
*  ===================
*
*  MAXITS  INTEGER, default = 5
*          The maximum number of iterations performed.
*
*  EXTRA   INTEGER, default = 2
*          The number of iterations performed after norm growth
*          criterion is satisfied, should be at least 1.
*
* =====================================================================
*
*     .. Parameters ..
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ),
     $                   CONE = ( 1.0E+0, 0.0E+0 ) )
      REAL               ZERO, ONE, TEN, ODM3, ODM1
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TEN = 1.0E+1,
     $                   ODM3 = 1.0E-3, ODM1 = 1.0E-1 )
      INTEGER            MAXITS, EXTRA
      PARAMETER          ( MAXITS = 5, EXTRA = 2 )
*     ..
*     .. Local Scalars ..
      INTEGER            B1, BLKSIZ, BN, GPIND, I, IINFO, INDRV1,
     $                   INDRV2, INDRV3, INDRV4, INDRV5, ITS, J, J1,
     $                   JBLK, JMAX, JR, NBLK, NRMCHK
      REAL               CTR, EPS, EPS1, NRM, ONENRM, ORTOL, PERTOL,
     $                   SCL, SEP, STPCRT, TOL, XJ, XJM
*     ..
*     .. Local Arrays ..
      INTEGER            ISEED( 4 )
*     ..
*     .. External Functions ..
      INTEGER            ISAMAX
      REAL               SASUM, SLAMCH, SNRM2
      EXTERNAL           ISAMAX, SASUM, SLAMCH, SNRM2
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SLAGTF, SLAGTS, SLARNV, SSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CMPLX, MAX, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      DO 10 I = 1, M
         IFAIL( I ) = 0
   10 CONTINUE
*
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 .OR. M.GT.N ) THEN
         INFO = -4
      ELSE IF( LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE
         DO 20 J = 2, M
            IF( IBLOCK( J ).LT.IBLOCK( J-1 ) ) THEN
               INFO = -6
               GO TO 30
            END IF
            IF( IBLOCK( J ).EQ.IBLOCK( J-1 ) .AND. W( J ).LT.W( J-1 ) )
     $           THEN
               INFO = -5
               GO TO 30
            END IF
   20    CONTINUE
   30    CONTINUE
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CSTEIN', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. M.EQ.0 ) THEN
         RETURN
      ELSE IF( N.EQ.1 ) THEN
         Z( 1, 1 ) = CONE
         RETURN
      END IF
*
*     Get machine constants.
*
      EPS = SLAMCH( 'Precision' )
*
*     Initialize seed for random number generator SLARNV.
*
      DO 40 I = 1, 4
         ISEED( I ) = 1
   40 CONTINUE
*
*     Initialize pointers.
*
      INDRV1 = 0
      INDRV2 = INDRV1 + N
      INDRV3 = INDRV2 + N
      INDRV4 = INDRV3 + N
      INDRV5 = INDRV4 + N
*
*     Compute eigenvectors of matrix blocks.
*
      J1 = 1
      DO 180 NBLK = 1, IBLOCK( M )
*
*        Find starting and ending indices of block nblk.
*
         IF( NBLK.EQ.1 ) THEN
            B1 = 1
         ELSE
            B1 = ISPLIT( NBLK-1 ) + 1
         END IF
         BN = ISPLIT( NBLK )
         BLKSIZ = BN - B1 + 1
         IF( BLKSIZ.EQ.1 )
     $      GO TO 60
         GPIND = B1
*
*        Compute reorthogonalization criterion and stopping criterion.
*
         ONENRM = ABS( D( B1 ) ) + ABS( E( B1 ) )
         ONENRM = MAX( ONENRM, ABS( D( BN ) )+ABS( E( BN-1 ) ) )
         DO 50 I = B1 + 1, BN - 1
            ONENRM = MAX( ONENRM, ABS( D( I ) )+ABS( E( I-1 ) )+
     $               ABS( E( I ) ) )
   50    CONTINUE
         ORTOL = ODM3*ONENRM
*
         STPCRT = SQRT( ODM1 / BLKSIZ )
*
*        Loop through eigenvalues of block nblk.
*
   60    CONTINUE
         JBLK = 0
         DO 170 J = J1, M
            IF( IBLOCK( J ).NE.NBLK ) THEN
               J1 = J
               GO TO 180
            END IF
            JBLK = JBLK + 1
            XJ = W( J )
*
*           Skip all the work if the block size is one.
*
            IF( BLKSIZ.EQ.1 ) THEN
               WORK( INDRV1+1 ) = ONE
               GO TO 140
            END IF
*
*           If eigenvalues j and j-1 are too close, add a relatively
*           small perturbation.
*
            IF( JBLK.GT.1 ) THEN
               EPS1 = ABS( EPS*XJ )
               PERTOL = TEN*EPS1
               SEP = XJ - XJM
               IF( SEP.LT.PERTOL )
     $            XJ = XJM + PERTOL
            END IF
*
            ITS = 0
            NRMCHK = 0
*
*           Get random starting vector.
*
            CALL SLARNV( 2, ISEED, BLKSIZ, WORK( INDRV1+1 ) )
*
*           Copy the matrix T so it won't be destroyed in factorization.
*
            CALL SCOPY( BLKSIZ, D( B1 ), 1, WORK( INDRV4+1 ), 1 )
            CALL SCOPY( BLKSIZ-1, E( B1 ), 1, WORK( INDRV2+2 ), 1 )
            CALL SCOPY( BLKSIZ-1, E( B1 ), 1, WORK( INDRV3+1 ), 1 )
*
*           Compute LU factors with partial pivoting  ( PT = LU )
*
            TOL = ZERO
            CALL SLAGTF( BLKSIZ, WORK( INDRV4+1 ), XJ, WORK( INDRV2+2 ),
     $                   WORK( INDRV3+1 ), TOL, WORK( INDRV5+1 ), IWORK,
     $                   IINFO )
*
*           Update iteration count.
*
   70       CONTINUE
            ITS = ITS + 1
            IF( ITS.GT.MAXITS )
     $         GO TO 120
*
*           Normalize and scale the righthand side vector Pb.
*
            SCL = BLKSIZ*ONENRM*MAX( EPS,
     $            ABS( WORK( INDRV4+BLKSIZ ) ) ) /
     $            SASUM( BLKSIZ, WORK( INDRV1+1 ), 1 )
            CALL SSCAL( BLKSIZ, SCL, WORK( INDRV1+1 ), 1 )
*
*           Solve the system LU = Pb.
*
            CALL SLAGTS( -1, BLKSIZ, WORK( INDRV4+1 ), WORK( INDRV2+2 ),
     $                   WORK( INDRV3+1 ), WORK( INDRV5+1 ), IWORK,
     $                   WORK( INDRV1+1 ), TOL, IINFO )
*
*           Reorthogonalize by modified Gram-Schmidt if eigenvalues are
*           close enough.
*
            IF( JBLK.EQ.1 )
     $         GO TO 110
            IF( ABS( XJ-XJM ).GT.ORTOL )
     $         GPIND = J
            IF( GPIND.NE.J ) THEN
               DO 100 I = GPIND, J - 1
                  CTR = ZERO
                  DO 80 JR = 1, BLKSIZ
                     CTR = CTR + WORK( INDRV1+JR )*
     $                     REAL( Z( B1-1+JR, I ) )
   80             CONTINUE
                  DO 90 JR = 1, BLKSIZ
                     WORK( INDRV1+JR ) = WORK( INDRV1+JR ) -
     $                                   CTR*REAL( Z( B1-1+JR, I ) )
   90             CONTINUE
  100          CONTINUE
            END IF
*
*           Check the infinity norm of the iterate.
*
  110       CONTINUE
            JMAX = ISAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
            NRM = ABS( WORK( INDRV1+JMAX ) )
*
*           Continue for additional iterations after norm reaches
*           stopping criterion.
*
            IF( NRM.LT.STPCRT )
     $         GO TO 70
            NRMCHK = NRMCHK + 1
            IF( NRMCHK.LT.EXTRA+1 )
     $         GO TO 70
*
            GO TO 130
*
*           If stopping criterion was not satisfied, update info and
*           store eigenvector number in array ifail.
*
  120       CONTINUE
            INFO = INFO + 1
            IFAIL( INFO ) = J
*
*           Accept iterate as jth eigenvector.
*
  130       CONTINUE
            SCL = ONE / SNRM2( BLKSIZ, WORK( INDRV1+1 ), 1 )
            JMAX = ISAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
            IF( WORK( INDRV1+JMAX ).LT.ZERO )
     $         SCL = -SCL
            CALL SSCAL( BLKSIZ, SCL, WORK( INDRV1+1 ), 1 )
  140       CONTINUE
            DO 150 I = 1, N
               Z( I, J ) = CZERO
  150       CONTINUE
            DO 160 I = 1, BLKSIZ
               Z( B1+I-1, J ) = CMPLX( WORK( INDRV1+I ), ZERO )
  160       CONTINUE
*
*           Save the shift to check eigenvalue spacing at next
*           iteration.
*
            XJM = XJ
*
  170    CONTINUE
  180 CONTINUE
*
      RETURN
*
*     End of CSTEIN
*
      END
