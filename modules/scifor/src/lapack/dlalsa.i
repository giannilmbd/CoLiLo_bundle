# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/dlalsa.f"
      SUBROUTINE DLALSA( ICOMPQ, SMLSIZ, N, NRHS, B, LDB, BX, LDBX, U,
     $                   LDU, VT, K, DIFL, DIFR, Z, POLES, GIVPTR,
     $                   GIVCOL, LDGCOL, PERM, GIVNUM, C, S, WORK,
     $                   IWORK, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            ICOMPQ, INFO, LDB, LDBX, LDGCOL, LDU, N, NRHS,
     $                   SMLSIZ
*     ..
*     .. Array Arguments ..
      INTEGER            GIVCOL( LDGCOL, * ), GIVPTR( * ), IWORK( * ),
     $                   K( * ), PERM( LDGCOL, * )
      DOUBLE PRECISION   B( LDB, * ), BX( LDBX, * ), C( * ),
     $                   DIFL( LDU, * ), DIFR( LDU, * ),
     $                   GIVNUM( LDU, * ), POLES( LDU, * ), S( * ),
     $                   U( LDU, * ), VT( LDU, * ), WORK( * ),
     $                   Z( LDU, * )
*     ..
*
*  Purpose
*  =======
*
*  DLALSA is an itermediate step in solving the least squares problem
*  by computing the SVD of the coefficient matrix in compact form (The
*  singular vectors are computed as products of simple orthorgonal
*  matrices.).
*
*  If ICOMPQ = 0, DLALSA applies the inverse of the left singular vector
*  matrix of an upper bidiagonal matrix to the right hand side; and if
*  ICOMPQ = 1, DLALSA applies the right singular vector matrix to the
*  right hand side. The singular vector matrices were generated in
*  compact form by DLALSA.
*
*  Arguments
*  =========
*
*
*  ICOMPQ (input) INTEGER
*         Specifies whether the left or the right singular vector
*         matrix is involved.
*         = 0: Left singular vector matrix
*         = 1: Right singular vector matrix
*
*  SMLSIZ (input) INTEGER
*         The maximum size of the subproblems at the bottom of the
*         computation tree.
*
*  N      (input) INTEGER
*         The row and column dimensions of the upper bidiagonal matrix.
*
*  NRHS   (input) INTEGER
*         The number of columns of B and BX. NRHS must be at least 1.
*
*  B      (input/output) DOUBLE PRECISION array, dimension ( LDB, NRHS )
*         On input, B contains the right hand sides of the least
*         squares problem in rows 1 through M.
*         On output, B contains the solution X in rows 1 through N.
*
*  LDB    (input) INTEGER
*         The leading dimension of B in the calling subprogram.
*         LDB must be at least max(1,MAX( M, N ) ).
*
*  BX     (output) DOUBLE PRECISION array, dimension ( LDBX, NRHS )
*         On exit, the result of applying the left or right singular
*         vector matrix to B.
*
*  LDBX   (input) INTEGER
*         The leading dimension of BX.
*
*  U      (input) DOUBLE PRECISION array, dimension ( LDU, SMLSIZ ).
*         On entry, U contains the left singular vector matrices of all
*         subproblems at the bottom level.
*
*  LDU    (input) INTEGER, LDU = > N.
*         The leading dimension of arrays U, VT, DIFL, DIFR,
*         POLES, GIVNUM, and Z.
*
*  VT     (input) DOUBLE PRECISION array, dimension ( LDU, SMLSIZ+1 ).
*         On entry, VT**T contains the right singular vector matrices of
*         all subproblems at the bottom level.
*
*  K      (input) INTEGER array, dimension ( N ).
*
*  DIFL   (input) DOUBLE PRECISION array, dimension ( LDU, NLVL ).
*         where NLVL = INT(log_2 (N/(SMLSIZ+1))) + 1.
*
*  DIFR   (input) DOUBLE PRECISION array, dimension ( LDU, 2 * NLVL ).
*         On entry, DIFL(*, I) and DIFR(*, 2 * I -1) record
*         distances between singular values on the I-th level and
*         singular values on the (I -1)-th level, and DIFR(*, 2 * I)
*         record the normalizing factors of the right singular vectors
*         matrices of subproblems on I-th level.
*
*  Z      (input) DOUBLE PRECISION array, dimension ( LDU, NLVL ).
*         On entry, Z(1, I) contains the components of the deflation-
*         adjusted updating row vector for subproblems on the I-th
*         level.
*
*  POLES  (input) DOUBLE PRECISION array, dimension ( LDU, 2 * NLVL ).
*         On entry, POLES(*, 2 * I -1: 2 * I) contains the new and old
*         singular values involved in the secular equations on the I-th
*         level.
*
*  GIVPTR (input) INTEGER array, dimension ( N ).
*         On entry, GIVPTR( I ) records the number of Givens
*         rotations performed on the I-th problem on the computation
*         tree.
*
*  GIVCOL (input) INTEGER array, dimension ( LDGCOL, 2 * NLVL ).
*         On entry, for each I, GIVCOL(*, 2 * I - 1: 2 * I) records the
*         locations of Givens rotations performed on the I-th level on
*         the computation tree.
*
*  LDGCOL (input) INTEGER, LDGCOL = > N.
*         The leading dimension of arrays GIVCOL and PERM.
*
*  PERM   (input) INTEGER array, dimension ( LDGCOL, NLVL ).
*         On entry, PERM(*, I) records permutations done on the I-th
*         level of the computation tree.
*
*  GIVNUM (input) DOUBLE PRECISION array, dimension ( LDU, 2 * NLVL ).
*         On entry, GIVNUM(*, 2 *I -1 : 2 * I) records the C- and S-
*         values of Givens rotations performed on the I-th level on the
*         computation tree.
*
*  C      (input) DOUBLE PRECISION array, dimension ( N ).
*         On entry, if the I-th subproblem is not square,
*         C( I ) contains the C-value of a Givens rotation related to
*         the right null space of the I-th subproblem.
*
*  S      (input) DOUBLE PRECISION array, dimension ( N ).
*         On entry, if the I-th subproblem is not square,
*         S( I ) contains the S-value of a Givens rotation related to
*         the right null space of the I-th subproblem.
*
*  WORK   (workspace) DOUBLE PRECISION array.
*         The dimension must be at least N.
*
*  IWORK  (workspace) INTEGER array.
*         The dimension must be at least 3 * N
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
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, I1, IC, IM1, INODE, J, LF, LL, LVL, LVL2,
     $                   ND, NDB1, NDIML, NDIMR, NL, NLF, NLP1, NLVL,
     $                   NR, NRF, NRP1, SQRE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DLALS0, DLASDT, XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) THEN
         INFO = -1
      ELSE IF( SMLSIZ.LT.3 ) THEN
         INFO = -2
      ELSE IF( N.LT.SMLSIZ ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.1 ) THEN
         INFO = -4
      ELSE IF( LDB.LT.N ) THEN
         INFO = -6
      ELSE IF( LDBX.LT.N ) THEN
         INFO = -8
      ELSE IF( LDU.LT.N ) THEN
         INFO = -10
      ELSE IF( LDGCOL.LT.N ) THEN
         INFO = -19
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLALSA', -INFO )
         RETURN
      END IF
*
*     Book-keeping and  setting up the computation tree.
*
      INODE = 1
      NDIML = INODE + N
      NDIMR = NDIML + N
*
      CALL DLASDT( N, NLVL, ND, IWORK( INODE ), IWORK( NDIML ),
     $             IWORK( NDIMR ), SMLSIZ )
*
*     The following code applies back the left singular vector factors.
*     For applying back the right singular vector factors, go to 50.
*
      IF( ICOMPQ.EQ.1 ) THEN
         GO TO 50
      END IF
*
*     The nodes on the bottom level of the tree were solved
*     by DLASDQ. The corresponding left and right singular vector
*     matrices are in explicit form. First apply back the left
*     singular vector matrices.
*
      NDB1 = ( ND+1 ) / 2
      DO 10 I = NDB1, ND
*
*        IC : center row of each node
*        NL : number of rows of left  subproblem
*        NR : number of rows of right subproblem
*        NLF: starting row of the left   subproblem
*        NRF: starting row of the right  subproblem
*
         I1 = I - 1
         IC = IWORK( INODE+I1 )
         NL = IWORK( NDIML+I1 )
         NR = IWORK( NDIMR+I1 )
         NLF = IC - NL
         NRF = IC + 1
         CALL DGEMM( 'T', 'N', NL, NRHS, NL, ONE, U( NLF, 1 ), LDU,
     $               B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )
         CALL DGEMM( 'T', 'N', NR, NRHS, NR, ONE, U( NRF, 1 ), LDU,
     $               B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )
   10 CONTINUE
*
*     Next copy the rows of B that correspond to unchanged rows
*     in the bidiagonal matrix to BX.
*
      DO 20 I = 1, ND
         IC = IWORK( INODE+I-1 )
         CALL DCOPY( NRHS, B( IC, 1 ), LDB, BX( IC, 1 ), LDBX )
   20 CONTINUE
*
*     Finally go through the left singular vector matrices of all
*     the other subproblems bottom-up on the tree.
*
      J = 2**NLVL
      SQRE = 0
*
      DO 40 LVL = NLVL, 1, -1
         LVL2 = 2*LVL - 1
*
*        find the first node LF and last node LL on
*        the current level LVL
*
         IF( LVL.EQ.1 ) THEN
            LF = 1
            LL = 1
         ELSE
            LF = 2**( LVL-1 )
            LL = 2*LF - 1
         END IF
         DO 30 I = LF, LL
            IM1 = I - 1
            IC = IWORK( INODE+IM1 )
            NL = IWORK( NDIML+IM1 )
            NR = IWORK( NDIMR+IM1 )
            NLF = IC - NL
            NRF = IC + 1
            J = J - 1
            CALL DLALS0( ICOMPQ, NL, NR, SQRE, NRHS, BX( NLF, 1 ), LDBX,
     $                   B( NLF, 1 ), LDB, PERM( NLF, LVL ),
     $                   GIVPTR( J ), GIVCOL( NLF, LVL2 ), LDGCOL,
     $                   GIVNUM( NLF, LVL2 ), LDU, POLES( NLF, LVL2 ),
     $                   DIFL( NLF, LVL ), DIFR( NLF, LVL2 ),
     $                   Z( NLF, LVL ), K( J ), C( J ), S( J ), WORK,
     $                   INFO )
   30    CONTINUE
   40 CONTINUE
      GO TO 90
*
*     ICOMPQ = 1: applying back the right singular vector factors.
*
   50 CONTINUE
*
*     First now go through the right singular vector matrices of all
*     the tree nodes top-down.
*
      J = 0
      DO 70 LVL = 1, NLVL
         LVL2 = 2*LVL - 1
*
*        Find the first node LF and last node LL on
*        the current level LVL.
*
         IF( LVL.EQ.1 ) THEN
            LF = 1
            LL = 1
         ELSE
            LF = 2**( LVL-1 )
            LL = 2*LF - 1
         END IF
         DO 60 I = LL, LF, -1
            IM1 = I - 1
            IC = IWORK( INODE+IM1 )
            NL = IWORK( NDIML+IM1 )
            NR = IWORK( NDIMR+IM1 )
            NLF = IC - NL
            NRF = IC + 1
            IF( I.EQ.LL ) THEN
               SQRE = 0
            ELSE
               SQRE = 1
            END IF
            J = J + 1
            CALL DLALS0( ICOMPQ, NL, NR, SQRE, NRHS, B( NLF, 1 ), LDB,
     $                   BX( NLF, 1 ), LDBX, PERM( NLF, LVL ),
     $                   GIVPTR( J ), GIVCOL( NLF, LVL2 ), LDGCOL,
     $                   GIVNUM( NLF, LVL2 ), LDU, POLES( NLF, LVL2 ),
     $                   DIFL( NLF, LVL ), DIFR( NLF, LVL2 ),
     $                   Z( NLF, LVL ), K( J ), C( J ), S( J ), WORK,
     $                   INFO )
   60    CONTINUE
   70 CONTINUE
*
*     The nodes on the bottom level of the tree were solved
*     by DLASDQ. The corresponding right singular vector
*     matrices are in explicit form. Apply them back.
*
      NDB1 = ( ND+1 ) / 2
      DO 80 I = NDB1, ND
         I1 = I - 1
         IC = IWORK( INODE+I1 )
         NL = IWORK( NDIML+I1 )
         NR = IWORK( NDIMR+I1 )
         NLP1 = NL + 1
         IF( I.EQ.ND ) THEN
            NRP1 = NR
         ELSE
            NRP1 = NR + 1
         END IF
         NLF = IC - NL
         NRF = IC + 1
         CALL DGEMM( 'T', 'N', NLP1, NRHS, NLP1, ONE, VT( NLF, 1 ), LDU,
     $               B( NLF, 1 ), LDB, ZERO, BX( NLF, 1 ), LDBX )
         CALL DGEMM( 'T', 'N', NRP1, NRHS, NRP1, ONE, VT( NRF, 1 ), LDU,
     $               B( NRF, 1 ), LDB, ZERO, BX( NRF, 1 ), LDBX )
   80 CONTINUE
*
   90 CONTINUE
*
      RETURN
*
*     End of DLALSA
*
      END
