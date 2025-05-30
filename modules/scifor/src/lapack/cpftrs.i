# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/cpftrs.f"
      SUBROUTINE CPFTRS( TRANSR, UPLO, N, NRHS, A, B, LDB, INFO )
*
*  -- LAPACK routine (version 3.3.1)                                    --
*
*  -- Contributed by Fred Gustavson of the IBM Watson Research Center --
*  -- April 2011                                                      --
*
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          TRANSR, UPLO
      INTEGER            INFO, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      COMPLEX            A( 0: * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  CPFTRS solves a system of linear equations A*X = B with a Hermitian
*  positive definite matrix A using the Cholesky factorization
*  A = U**H*U or A = L*L**H computed by CPFTRF.
*
*  Arguments
*  =========
*
*  TRANSR    (input) CHARACTER*1
*          = 'N':  The Normal TRANSR of RFP A is stored;
*          = 'C':  The Conjugate-transpose TRANSR of RFP A is stored.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of RFP A is stored;
*          = 'L':  Lower triangle of RFP A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) COMPLEX array, dimension ( N*(N+1)/2 );
*          The triangular factor U or L from the Cholesky factorization
*          of RFP A = U**H*U or RFP A = L*L**H, as computed by CPFTRF.
*          See note below for more details about RFP A.
*
*  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  We first consider Standard Packed Format when N is even.
*  We give an example where N = 6.
*
*      AP is Upper             AP is Lower
*
*   00 01 02 03 04 05       00
*      11 12 13 14 15       10 11
*         22 23 24 25       20 21 22
*            33 34 35       30 31 32 33
*               44 45       40 41 42 43 44
*                  55       50 51 52 53 54 55
*
*
*  Let TRANSR = 'N'. RFP holds AP as follows:
*  For UPLO = 'U' the upper trapezoid A(0:5,0:2) consists of the last
*  three columns of AP upper. The lower triangle A(4:6,0:2) consists of
*  conjugate-transpose of the first three columns of AP upper.
*  For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first
*  three columns of AP lower. The upper triangle A(0:2,0:2) consists of
*  conjugate-transpose of the last three columns of AP lower.
*  To denote conjugate we place -- above the element. This covers the
*  case N even and TRANSR = 'N'.
*
*         RFP A                   RFP A
*
*                                -- -- --
*        03 04 05                33 43 53
*                                   -- --
*        13 14 15                00 44 54
*                                      --
*        23 24 25                10 11 55
*
*        33 34 35                20 21 22
*        --
*        00 44 45                30 31 32
*        -- --
*        01 11 55                40 41 42
*        -- -- --
*        02 12 22                50 51 52
*
*  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate-
*  transpose of RFP A above. One therefore gets:
*
*
*           RFP A                   RFP A
*
*     -- -- -- --                -- -- -- -- -- --
*     03 13 23 33 00 01 02    33 00 10 20 30 40 50
*     -- -- -- -- --                -- -- -- -- --
*     04 14 24 34 44 11 12    43 44 11 21 31 41 51
*     -- -- -- -- -- --                -- -- -- --
*     05 15 25 35 45 55 22    53 54 55 22 32 42 52
*
*
*  We next  consider Standard Packed Format when N is odd.
*  We give an example where N = 5.
*
*     AP is Upper                 AP is Lower
*
*   00 01 02 03 04              00
*      11 12 13 14              10 11
*         22 23 24              20 21 22
*            33 34              30 31 32 33
*               44              40 41 42 43 44
*
*
*  Let TRANSR = 'N'. RFP holds AP as follows:
*  For UPLO = 'U' the upper trapezoid A(0:4,0:2) consists of the last
*  three columns of AP upper. The lower triangle A(3:4,0:1) consists of
*  conjugate-transpose of the first two   columns of AP upper.
*  For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first
*  three columns of AP lower. The upper triangle A(0:1,1:2) consists of
*  conjugate-transpose of the last two   columns of AP lower.
*  To denote conjugate we place -- above the element. This covers the
*  case N odd  and TRANSR = 'N'.
*
*         RFP A                   RFP A
*
*                                   -- --
*        02 03 04                00 33 43
*                                      --
*        12 13 14                10 11 44
*
*        22 23 24                20 21 22
*        --
*        00 33 34                30 31 32
*        -- --
*        01 11 44                40 41 42
*
*  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate-
*  transpose of RFP A above. One therefore gets:
*
*
*           RFP A                   RFP A
*
*     -- -- --                   -- -- -- -- -- --
*     02 12 22 00 01             00 10 20 30 40 50
*     -- -- -- --                   -- -- -- -- --
*     03 13 23 33 11             33 11 21 31 41 51
*     -- -- -- -- --                   -- -- -- --
*     04 14 24 34 44             43 44 22 32 42 52
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, NORMALTRANSR
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, CTFSM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NORMALTRANSR = LSAME( TRANSR, 'N' )
      LOWER = LSAME( UPLO, 'L' )
      IF( .NOT.NORMALTRANSR .AND. .NOT.LSAME( TRANSR, 'C' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CPFTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
*     start execution: there are two triangular solves
*
      IF( LOWER ) THEN
         CALL CTFSM( TRANSR, 'L', UPLO, 'N', 'N', N, NRHS, CONE, A, B,
     $               LDB )
         CALL CTFSM( TRANSR, 'L', UPLO, 'C', 'N', N, NRHS, CONE, A, B,
     $               LDB )
      ELSE
         CALL CTFSM( TRANSR, 'L', UPLO, 'C', 'N', N, NRHS, CONE, A, B,
     $               LDB )
         CALL CTFSM( TRANSR, 'L', UPLO, 'N', 'N', N, NRHS, CONE, A, B,
     $               LDB )
      END IF
*
      RETURN
*
*     End of CPFTRS
*
      END
