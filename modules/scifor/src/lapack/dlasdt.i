# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/dlasdt.f"
      SUBROUTINE DLASDT( N, LVL, ND, INODE, NDIML, NDIMR, MSUB )
*
*  -- LAPACK auxiliary routine (version 3.2.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2010
*
*     .. Scalar Arguments ..
      INTEGER            LVL, MSUB, N, ND
*     ..
*     .. Array Arguments ..
      INTEGER            INODE( * ), NDIML( * ), NDIMR( * )
*     ..
*
*  Purpose
*  =======
*
*  DLASDT creates a tree of subproblems for bidiagonal divide and
*  conquer.
*
*  Arguments
*  =========
*
*   N      (input) INTEGER
*          On entry, the number of diagonal elements of the
*          bidiagonal matrix.
*
*   LVL    (output) INTEGER
*          On exit, the number of levels on the computation tree.
*
*   ND     (output) INTEGER
*          On exit, the number of nodes on the tree.
*
*   INODE  (output) INTEGER array, dimension ( N )
*          On exit, centers of subproblems.
*
*   NDIML  (output) INTEGER array, dimension ( N )
*          On exit, row dimensions of left children.
*
*   NDIMR  (output) INTEGER array, dimension ( N )
*          On exit, row dimensions of right children.
*
*   MSUB   (input) INTEGER
*          On entry, the maximum row dimension each subproblem at the
*          bottom of the tree can be of.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Ming Gu and Huan Ren, Computer Science Division, University of
*     California at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IL, IR, LLST, MAXN, NCRNT, NLVL
      DOUBLE PRECISION   TEMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, LOG, MAX
*     ..
*     .. Executable Statements ..
*
*     Find the number of levels on the tree.
*
      MAXN = MAX( 1, N )
      TEMP = LOG( DBLE( MAXN ) / DBLE( MSUB+1 ) ) / LOG( TWO )
      LVL = INT( TEMP ) + 1
*
      I = N / 2
      INODE( 1 ) = I + 1
      NDIML( 1 ) = I
      NDIMR( 1 ) = N - I - 1
      IL = 0
      IR = 1
      LLST = 1
      DO 20 NLVL = 1, LVL - 1
*
*        Constructing the tree at (NLVL+1)-st level. The number of
*        nodes created on this level is LLST * 2.
*
         DO 10 I = 0, LLST - 1
            IL = IL + 2
            IR = IR + 2
            NCRNT = LLST + I
            NDIML( IL ) = NDIML( NCRNT ) / 2
            NDIMR( IL ) = NDIML( NCRNT ) - NDIML( IL ) - 1
            INODE( IL ) = INODE( NCRNT ) - NDIMR( IL ) - 1
            NDIML( IR ) = NDIMR( NCRNT ) / 2
            NDIMR( IR ) = NDIMR( NCRNT ) - NDIML( IR ) - 1
            INODE( IR ) = INODE( NCRNT ) + NDIML( IR ) + 1
   10    CONTINUE
         LLST = LLST*2
   20 CONTINUE
      ND = LLST*2 - 1
*
      RETURN
*
*     End of DLASDT
*
      END
