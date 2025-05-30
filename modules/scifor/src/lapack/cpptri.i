# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/cpptri.f"
      SUBROUTINE CPPTRI( UPLO, N, AP, INFO )
*
*  -- LAPACK routine (version 3.3.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- April 2011                                                      --
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      COMPLEX            AP( * )
*     ..
*
*  Purpose
*  =======
*
*  CPPTRI computes the inverse of a complex Hermitian positive definite
*  matrix A using the Cholesky factorization A = U**H*U or A = L*L**H
*  computed by CPPTRF.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangular factor is stored in AP;
*          = 'L':  Lower triangular factor is stored in AP.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  AP      (input/output) COMPLEX array, dimension (N*(N+1)/2)
*          On entry, the triangular factor U or L from the Cholesky
*          factorization A = U**H*U or A = L*L**H, packed columnwise as
*          a linear array.  The j-th column of U or L is stored in the
*          array AP as follows:
*          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;
*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.
*
*          On exit, the upper or lower triangle of the (Hermitian)
*          inverse of A, overwriting the input factor U or L.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the (i,i) element of the factor U or L is
*                zero, and the inverse could not be computed.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, JC, JJ, JJN
      REAL               AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX            CDOTC
      EXTERNAL           LSAME, CDOTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHPR, CSSCAL, CTPMV, CTPTRI, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CPPTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Invert the triangular Cholesky factor U or L.
*
      CALL CTPTRI( UPLO, 'Non-unit', N, AP, INFO )
      IF( INFO.GT.0 )
     $   RETURN
      IF( UPPER ) THEN
*
*        Compute the product inv(U) * inv(U)**H.
*
         JJ = 0
         DO 10 J = 1, N
            JC = JJ + 1
            JJ = JJ + J
            IF( J.GT.1 )
     $         CALL CHPR( 'Upper', J-1, ONE, AP( JC ), 1, AP )
            AJJ = AP( JJ )
            CALL CSSCAL( J, AJJ, AP( JC ), 1 )
   10    CONTINUE
*
      ELSE
*
*        Compute the product inv(L)**H * inv(L).
*
         JJ = 1
         DO 20 J = 1, N
            JJN = JJ + N - J + 1
            AP( JJ ) = REAL( CDOTC( N-J+1, AP( JJ ), 1, AP( JJ ), 1 ) )
            IF( J.LT.N )
     $         CALL CTPMV( 'Lower', 'Conjugate transpose', 'Non-unit',
     $                     N-J, AP( JJN ), AP( JJ+1 ), 1 )
            JJ = JJN
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of CPPTRI
*
      END
