# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/ctgexc.f"
      SUBROUTINE CTGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,
     $                   LDZ, IFST, ILST, INFO )
*
*  -- LAPACK routine (version 3.3.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- April 2011                                                      --
*
*     .. Scalar Arguments ..
      LOGICAL            WANTQ, WANTZ
      INTEGER            IFST, ILST, INFO, LDA, LDB, LDQ, LDZ, N
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
     $                   Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  CTGEXC reorders the generalized Schur decomposition of a complex
*  matrix pair (A,B), using an unitary equivalence transformation
*  (A, B) := Q * (A, B) * Z**H, so that the diagonal block of (A, B) with
*  row index IFST is moved to row ILST.
*
*  (A, B) must be in generalized Schur canonical form, that is, A and
*  B are both upper triangular.
*
*  Optionally, the matrices Q and Z of generalized Schur vectors are
*  updated.
*
*         Q(in) * A(in) * Z(in)**H = Q(out) * A(out) * Z(out)**H
*         Q(in) * B(in) * Z(in)**H = Q(out) * B(out) * Z(out)**H
*
*  Arguments
*  =========
*
*  WANTQ   (input) LOGICAL
*          .TRUE. : update the left transformation matrix Q;
*          .FALSE.: do not update Q.
*
*  WANTZ   (input) LOGICAL
*          .TRUE. : update the right transformation matrix Z;
*          .FALSE.: do not update Z.
*
*  N       (input) INTEGER
*          The order of the matrices A and B. N >= 0.
*
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the upper triangular matrix A in the pair (A, B).
*          On exit, the updated matrix A.
*
*  LDA     (input)  INTEGER
*          The leading dimension of the array A. LDA >= max(1,N).
*
*  B       (input/output) COMPLEX array, dimension (LDB,N)
*          On entry, the upper triangular matrix B in the pair (A, B).
*          On exit, the updated matrix B.
*
*  LDB     (input)  INTEGER
*          The leading dimension of the array B. LDB >= max(1,N).
*
*  Q       (input/output) COMPLEX array, dimension (LDZ,N)
*          On entry, if WANTQ = .TRUE., the unitary matrix Q.
*          On exit, the updated matrix Q.
*          If WANTQ = .FALSE., Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q. LDQ >= 1;
*          If WANTQ = .TRUE., LDQ >= N.
*
*  Z       (input/output) COMPLEX array, dimension (LDZ,N)
*          On entry, if WANTZ = .TRUE., the unitary matrix Z.
*          On exit, the updated matrix Z.
*          If WANTZ = .FALSE., Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z. LDZ >= 1;
*          If WANTZ = .TRUE., LDZ >= N.
*
*  IFST    (input) INTEGER
*  ILST    (input/output) INTEGER
*          Specify the reordering of the diagonal blocks of (A, B).
*          The block with row index IFST is moved to row ILST, by a
*          sequence of swapping between adjacent blocks.
*
*  INFO    (output) INTEGER
*           =0:  Successful exit.
*           <0:  if INFO = -i, the i-th argument had an illegal value.
*           =1:  The transformed matrix pair (A, B) would be too far
*                from generalized Schur form; the problem is ill-
*                conditioned. (A, B) may have been partially reordered,
*                and ILST points to the first row of the current
*                position of the block being moved.
*
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
*     Umea University, S-901 87 Umea, Sweden.
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
*      Sweden, December 1993, Revised April 1994, Also as LAPACK working
*      Note 75. To appear in ACM Trans. on Math. Software, Vol 22, No 1,
*      1996.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            HERE
*     ..
*     .. External Subroutines ..
      EXTERNAL           CTGEX2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Decode and test input arguments.
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDQ.LT.1 .OR. WANTQ .AND. ( LDQ.LT.MAX( 1, N ) ) ) THEN
         INFO = -9
      ELSE IF( LDZ.LT.1 .OR. WANTZ .AND. ( LDZ.LT.MAX( 1, N ) ) ) THEN
         INFO = -11
      ELSE IF( IFST.LT.1 .OR. IFST.GT.N ) THEN
         INFO = -12
      ELSE IF( ILST.LT.1 .OR. ILST.GT.N ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTGEXC', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
      IF( IFST.EQ.ILST )
     $   RETURN
*
      IF( IFST.LT.ILST ) THEN
*
         HERE = IFST
*
   10    CONTINUE
*
*        Swap with next one below
*
         CALL CTGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ,
     $                HERE, INFO )
         IF( INFO.NE.0 ) THEN
            ILST = HERE
            RETURN
         END IF
         HERE = HERE + 1
         IF( HERE.LT.ILST )
     $      GO TO 10
         HERE = HERE - 1
      ELSE
         HERE = IFST - 1
*
   20    CONTINUE
*
*        Swap with next one above
*
         CALL CTGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ,
     $                HERE, INFO )
         IF( INFO.NE.0 ) THEN
            ILST = HERE
            RETURN
         END IF
         HERE = HERE - 1
         IF( HERE.GE.ILST )
     $      GO TO 20
         HERE = HERE + 1
      END IF
      ILST = HERE
      RETURN
*
*     End of CTGEXC
*
      END
