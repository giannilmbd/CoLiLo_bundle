# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/dlals0.f"
      SUBROUTINE DLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B, LDB, BX, LDBX,
     $                   PERM, GIVPTR, GIVCOL, LDGCOL, GIVNUM, LDGNUM,
     $                   POLES, DIFL, DIFR, Z, K, C, S, WORK, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            GIVPTR, ICOMPQ, INFO, K, LDB, LDBX, LDGCOL,
     $                   LDGNUM, NL, NR, NRHS, SQRE
      DOUBLE PRECISION   C, S
*     ..
*     .. Array Arguments ..
      INTEGER            GIVCOL( LDGCOL, * ), PERM( * )
      DOUBLE PRECISION   B( LDB, * ), BX( LDBX, * ), DIFL( * ),
     $                   DIFR( LDGNUM, * ), GIVNUM( LDGNUM, * ),
     $                   POLES( LDGNUM, * ), WORK( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  DLALS0 applies back the multiplying factors of either the left or the
*  right singular vector matrix of a diagonal matrix appended by a row
*  to the right hand side matrix B in solving the least squares problem
*  using the divide-and-conquer SVD approach.
*
*  For the left singular vector matrix, three types of orthogonal
*  matrices are involved:
*
*  (1L) Givens rotations: the number of such rotations is GIVPTR; the
*       pairs of columns/rows they were applied to are stored in GIVCOL;
*       and the C- and S-values of these rotations are stored in GIVNUM.
*
*  (2L) Permutation. The (NL+1)-st row of B is to be moved to the first
*       row, and for J=2:N, PERM(J)-th row of B is to be moved to the
*       J-th row.
*
*  (3L) The left singular vector matrix of the remaining matrix.
*
*  For the right singular vector matrix, four types of orthogonal
*  matrices are involved:
*
*  (1R) The right singular vector matrix of the remaining matrix.
*
*  (2R) If SQRE = 1, one extra Givens rotation to generate the right
*       null space.
*
*  (3R) The inverse transformation of (2L).
*
*  (4R) The inverse transformation of (1L).
*
*  Arguments
*  =========
*
*  ICOMPQ (input) INTEGER
*         Specifies whether singular vectors are to be computed in
*         factored form:
*         = 0: Left singular vector matrix.
*         = 1: Right singular vector matrix.
*
*  NL     (input) INTEGER
*         The row dimension of the upper block. NL >= 1.
*
*  NR     (input) INTEGER
*         The row dimension of the lower block. NR >= 1.
*
*  SQRE   (input) INTEGER
*         = 0: the lower block is an NR-by-NR square matrix.
*         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.
*
*         The bidiagonal matrix has row dimension N = NL + NR + 1,
*         and column dimension M = N + SQRE.
*
*  NRHS   (input) INTEGER
*         The number of columns of B and BX. NRHS must be at least 1.
*
*  B      (input/output) DOUBLE PRECISION array, dimension ( LDB, NRHS )
*         On input, B contains the right hand sides of the least
*         squares problem in rows 1 through M. On output, B contains
*         the solution X in rows 1 through N.
*
*  LDB    (input) INTEGER
*         The leading dimension of B. LDB must be at least
*         max(1,MAX( M, N ) ).
*
*  BX     (workspace) DOUBLE PRECISION array, dimension ( LDBX, NRHS )
*
*  LDBX   (input) INTEGER
*         The leading dimension of BX.
*
*  PERM   (input) INTEGER array, dimension ( N )
*         The permutations (from deflation and sorting) applied
*         to the two blocks.
*
*  GIVPTR (input) INTEGER
*         The number of Givens rotations which took place in this
*         subproblem.
*
*  GIVCOL (input) INTEGER array, dimension ( LDGCOL, 2 )
*         Each pair of numbers indicates a pair of rows/columns
*         involved in a Givens rotation.
*
*  LDGCOL (input) INTEGER
*         The leading dimension of GIVCOL, must be at least N.
*
*  GIVNUM (input) DOUBLE PRECISION array, dimension ( LDGNUM, 2 )
*         Each number indicates the C or S value used in the
*         corresponding Givens rotation.
*
*  LDGNUM (input) INTEGER
*         The leading dimension of arrays DIFR, POLES and
*         GIVNUM, must be at least K.
*
*  POLES  (input) DOUBLE PRECISION array, dimension ( LDGNUM, 2 )
*         On entry, POLES(1:K, 1) contains the new singular
*         values obtained from solving the secular equation, and
*         POLES(1:K, 2) is an array containing the poles in the secular
*         equation.
*
*  DIFL   (input) DOUBLE PRECISION array, dimension ( K ).
*         On entry, DIFL(I) is the distance between I-th updated
*         (undeflated) singular value and the I-th (undeflated) old
*         singular value.
*
*  DIFR   (input) DOUBLE PRECISION array, dimension ( LDGNUM, 2 ).
*         On entry, DIFR(I, 1) contains the distances between I-th
*         updated (undeflated) singular value and the I+1-th
*         (undeflated) old singular value. And DIFR(I, 2) is the
*         normalizing factor for the I-th right singular vector.
*
*  Z      (input) DOUBLE PRECISION array, dimension ( K )
*         Contain the components of the deflation-adjusted updating row
*         vector.
*
*  K      (input) INTEGER
*         Contains the dimension of the non-deflated matrix,
*         This is the order of the related secular equation. 1 <= K <=N.
*
*  C      (input) DOUBLE PRECISION
*         C contains garbage if SQRE =0 and the C-value of a Givens
*         rotation related to the right null space if SQRE = 1.
*
*  S      (input) DOUBLE PRECISION
*         S contains garbage if SQRE =0 and the S-value of a Givens
*         rotation related to the right null space if SQRE = 1.
*
*  WORK   (workspace) DOUBLE PRECISION array, dimension ( K )
*
*  INFO   (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Ming Gu and Ren-Cang Li, Computer Science Division, University of
*       California at Berkeley, USA
*     Osni Marques, LBNL/NERSC, USA
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO, NEGONE
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0, NEGONE = -1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, M, N, NLP1
      DOUBLE PRECISION   DIFLJ, DIFRJ, DJ, DSIGJ, DSIGJP, TEMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMV, DLACPY, DLASCL, DROT, DSCAL,
     $                   XERBLA
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3, DNRM2
      EXTERNAL           DLAMC3, DNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) THEN
         INFO = -1
      ELSE IF( NL.LT.1 ) THEN
         INFO = -2
      ELSE IF( NR.LT.1 ) THEN
         INFO = -3
      ELSE IF( ( SQRE.LT.0 ) .OR. ( SQRE.GT.1 ) ) THEN
         INFO = -4
      END IF
*
      N = NL + NR + 1
*
      IF( NRHS.LT.1 ) THEN
         INFO = -5
      ELSE IF( LDB.LT.N ) THEN
         INFO = -7
      ELSE IF( LDBX.LT.N ) THEN
         INFO = -9
      ELSE IF( GIVPTR.LT.0 ) THEN
         INFO = -11
      ELSE IF( LDGCOL.LT.N ) THEN
         INFO = -13
      ELSE IF( LDGNUM.LT.N ) THEN
         INFO = -15
      ELSE IF( K.LT.1 ) THEN
         INFO = -20
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLALS0', -INFO )
         RETURN
      END IF
*
      M = N + SQRE
      NLP1 = NL + 1
*
      IF( ICOMPQ.EQ.0 ) THEN
*
*        Apply back orthogonal transformations from the left.
*
*        Step (1L): apply back the Givens rotations performed.
*
         DO 10 I = 1, GIVPTR
            CALL DROT( NRHS, B( GIVCOL( I, 2 ), 1 ), LDB,
     $                 B( GIVCOL( I, 1 ), 1 ), LDB, GIVNUM( I, 2 ),
     $                 GIVNUM( I, 1 ) )
   10    CONTINUE
*
*        Step (2L): permute rows of B.
*
         CALL DCOPY( NRHS, B( NLP1, 1 ), LDB, BX( 1, 1 ), LDBX )
         DO 20 I = 2, N
            CALL DCOPY( NRHS, B( PERM( I ), 1 ), LDB, BX( I, 1 ), LDBX )
   20    CONTINUE
*
*        Step (3L): apply the inverse of the left singular vector
*        matrix to BX.
*
         IF( K.EQ.1 ) THEN
            CALL DCOPY( NRHS, BX, LDBX, B, LDB )
            IF( Z( 1 ).LT.ZERO ) THEN
               CALL DSCAL( NRHS, NEGONE, B, LDB )
            END IF
         ELSE
            DO 50 J = 1, K
               DIFLJ = DIFL( J )
               DJ = POLES( J, 1 )
               DSIGJ = -POLES( J, 2 )
               IF( J.LT.K ) THEN
                  DIFRJ = -DIFR( J, 1 )
                  DSIGJP = -POLES( J+1, 2 )
               END IF
               IF( ( Z( J ).EQ.ZERO ) .OR. ( POLES( J, 2 ).EQ.ZERO ) )
     $              THEN
                  WORK( J ) = ZERO
               ELSE
                  WORK( J ) = -POLES( J, 2 )*Z( J ) / DIFLJ /
     $                        ( POLES( J, 2 )+DJ )
               END IF
               DO 30 I = 1, J - 1
                  IF( ( Z( I ).EQ.ZERO ) .OR.
     $                ( POLES( I, 2 ).EQ.ZERO ) ) THEN
                     WORK( I ) = ZERO
                  ELSE
                     WORK( I ) = POLES( I, 2 )*Z( I ) /
     $                           ( DLAMC3( POLES( I, 2 ), DSIGJ )-
     $                           DIFLJ ) / ( POLES( I, 2 )+DJ )
                  END IF
   30          CONTINUE
               DO 40 I = J + 1, K
                  IF( ( Z( I ).EQ.ZERO ) .OR.
     $                ( POLES( I, 2 ).EQ.ZERO ) ) THEN
                     WORK( I ) = ZERO
                  ELSE
                     WORK( I ) = POLES( I, 2 )*Z( I ) /
     $                           ( DLAMC3( POLES( I, 2 ), DSIGJP )+
     $                           DIFRJ ) / ( POLES( I, 2 )+DJ )
                  END IF
   40          CONTINUE
               WORK( 1 ) = NEGONE
               TEMP = DNRM2( K, WORK, 1 )
               CALL DGEMV( 'T', K, NRHS, ONE, BX, LDBX, WORK, 1, ZERO,
     $                     B( J, 1 ), LDB )
               CALL DLASCL( 'G', 0, 0, TEMP, ONE, 1, NRHS, B( J, 1 ),
     $                      LDB, INFO )
   50       CONTINUE
         END IF
*
*        Move the deflated rows of BX to B also.
*
         IF( K.LT.MAX( M, N ) )
     $      CALL DLACPY( 'A', N-K, NRHS, BX( K+1, 1 ), LDBX,
     $                   B( K+1, 1 ), LDB )
      ELSE
*
*        Apply back the right orthogonal transformations.
*
*        Step (1R): apply back the new right singular vector matrix
*        to B.
*
         IF( K.EQ.1 ) THEN
            CALL DCOPY( NRHS, B, LDB, BX, LDBX )
         ELSE
            DO 80 J = 1, K
               DSIGJ = POLES( J, 2 )
               IF( Z( J ).EQ.ZERO ) THEN
                  WORK( J ) = ZERO
               ELSE
                  WORK( J ) = -Z( J ) / DIFL( J ) /
     $                        ( DSIGJ+POLES( J, 1 ) ) / DIFR( J, 2 )
               END IF
               DO 60 I = 1, J - 1
                  IF( Z( J ).EQ.ZERO ) THEN
                     WORK( I ) = ZERO
                  ELSE
                     WORK( I ) = Z( J ) / ( DLAMC3( DSIGJ, -POLES( I+1,
     $                           2 ) )-DIFR( I, 1 ) ) /
     $                           ( DSIGJ+POLES( I, 1 ) ) / DIFR( I, 2 )
                  END IF
   60          CONTINUE
               DO 70 I = J + 1, K
                  IF( Z( J ).EQ.ZERO ) THEN
                     WORK( I ) = ZERO
                  ELSE
                     WORK( I ) = Z( J ) / ( DLAMC3( DSIGJ, -POLES( I,
     $                           2 ) )-DIFL( I ) ) /
     $                           ( DSIGJ+POLES( I, 1 ) ) / DIFR( I, 2 )
                  END IF
   70          CONTINUE
               CALL DGEMV( 'T', K, NRHS, ONE, B, LDB, WORK, 1, ZERO,
     $                     BX( J, 1 ), LDBX )
   80       CONTINUE
         END IF
*
*        Step (2R): if SQRE = 1, apply back the rotation that is
*        related to the right null space of the subproblem.
*
         IF( SQRE.EQ.1 ) THEN
            CALL DCOPY( NRHS, B( M, 1 ), LDB, BX( M, 1 ), LDBX )
            CALL DROT( NRHS, BX( 1, 1 ), LDBX, BX( M, 1 ), LDBX, C, S )
         END IF
         IF( K.LT.MAX( M, N ) )
     $      CALL DLACPY( 'A', N-K, NRHS, B( K+1, 1 ), LDB, BX( K+1, 1 ),
     $                   LDBX )
*
*        Step (3R): permute rows of B.
*
         CALL DCOPY( NRHS, BX( 1, 1 ), LDBX, B( NLP1, 1 ), LDB )
         IF( SQRE.EQ.1 ) THEN
            CALL DCOPY( NRHS, BX( M, 1 ), LDBX, B( M, 1 ), LDB )
         END IF
         DO 90 I = 2, N
            CALL DCOPY( NRHS, BX( I, 1 ), LDBX, B( PERM( I ), 1 ), LDB )
   90    CONTINUE
*
*        Step (4R): apply back the Givens rotations performed.
*
         DO 100 I = GIVPTR, 1, -1
            CALL DROT( NRHS, B( GIVCOL( I, 2 ), 1 ), LDB,
     $                 B( GIVCOL( I, 1 ), 1 ), LDB, GIVNUM( I, 2 ),
     $                 -GIVNUM( I, 1 ) )
  100    CONTINUE
      END IF
*
      RETURN
*
*     End of DLALS0
*
      END
