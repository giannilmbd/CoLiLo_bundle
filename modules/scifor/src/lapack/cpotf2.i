# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/cpotf2.f"
      SUBROUTINE CPOTF2( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 3.3.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- April 2011                                                      --
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  CPOTF2 computes the Cholesky factorization of a complex Hermitian
*  positive definite matrix A.
*
*  The factorization has the form
*     A = U**H * U ,  if UPLO = 'U', or
*     A = L  * L**H,  if UPLO = 'L',
*  where U is an upper triangular matrix and L is lower triangular.
*
*  This is the unblocked version of the algorithm, calling Level 2 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          Hermitian matrix A is stored.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
*          n by n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n by n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if INFO = 0, the factor U or L from the Cholesky
*          factorization A = U**H *U  or A = L*L**H.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, the leading minor of order k is not
*               positive definite, and the factorization could not be
*               completed.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J
      REAL               AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, SISNAN
      COMPLEX            CDOTC
      EXTERNAL           LSAME, CDOTC, SISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMV, CLACGV, CSSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, REAL, SQRT
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
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CPOTF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Compute the Cholesky factorization A = U**H *U.
*
         DO 10 J = 1, N
*
*           Compute U(J,J) and test for non-positive-definiteness.
*
            AJJ = REAL( A( J, J ) ) - CDOTC( J-1, A( 1, J ), 1,
     $            A( 1, J ), 1 )
            IF( AJJ.LE.ZERO.OR.SISNAN( AJJ ) ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
*
*           Compute elements J+1:N of row J.
*
            IF( J.LT.N ) THEN
               CALL CLACGV( J-1, A( 1, J ), 1 )
               CALL CGEMV( 'Transpose', J-1, N-J, -CONE, A( 1, J+1 ),
     $                     LDA, A( 1, J ), 1, CONE, A( J, J+1 ), LDA )
               CALL CLACGV( J-1, A( 1, J ), 1 )
               CALL CSSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
            END IF
   10    CONTINUE
      ELSE
*
*        Compute the Cholesky factorization A = L*L**H.
*
         DO 20 J = 1, N
*
*           Compute L(J,J) and test for non-positive-definiteness.
*
            AJJ = REAL( A( J, J ) ) - CDOTC( J-1, A( J, 1 ), LDA,
     $            A( J, 1 ), LDA )
            IF( AJJ.LE.ZERO.OR.SISNAN( AJJ ) ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
*
*           Compute elements J+1:N of column J.
*
            IF( J.LT.N ) THEN
               CALL CLACGV( J-1, A( J, 1 ), LDA )
               CALL CGEMV( 'No transpose', N-J, J-1, -CONE, A( J+1, 1 ),
     $                     LDA, A( J, 1 ), LDA, CONE, A( J+1, J ), 1 )
               CALL CLACGV( J-1, A( J, 1 ), LDA )
               CALL CSSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
      GO TO 40
*
   30 CONTINUE
      INFO = J
*
   40 CONTINUE
      RETURN
*
*     End of CPOTF2
*
      END
