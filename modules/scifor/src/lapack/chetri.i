# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/chetri.f"
      SUBROUTINE CHETRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
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
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CHETRI computes the inverse of a complex Hermitian indefinite matrix
*  A using the factorization A = U*D*U**H or A = L*D*L**H computed by
*  CHETRF.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the details of the factorization are stored
*          as an upper or lower triangular matrix.
*          = 'U':  Upper triangular, form is A = U*D*U**H;
*          = 'L':  Lower triangular, form is A = L*D*L**H.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the block diagonal matrix D and the multipliers
*          used to obtain the factor U or L as computed by CHETRF.
*
*          On exit, if INFO = 0, the (Hermitian) inverse of the original
*          matrix.  If UPLO = 'U', the upper triangular part of the
*          inverse is formed and the part of A below the diagonal is not
*          referenced; if UPLO = 'L' the lower triangular part of the
*          inverse is formed and the part of A above the diagonal is
*          not referenced.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          Details of the interchanges and the block structure of D
*          as determined by CHETRF.
*
*  WORK    (workspace) COMPLEX array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
*               inverse could not be computed.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      COMPLEX            CONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, CONE = ( 1.0E+0, 0.0E+0 ),
     $                   ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, K, KP, KSTEP
      REAL               AK, AKP1, D, T
      COMPLEX            AKKP1, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX            CDOTC
      EXTERNAL           LSAME, CDOTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           CCOPY, CHEMV, CSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CONJG, MAX, REAL
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
         CALL XERBLA( 'CHETRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Check that the diagonal matrix D is nonsingular.
*
      IF( UPPER ) THEN
*
*        Upper triangular storage: examine D from bottom to top
*
         DO 10 INFO = N, 1, -1
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   10    CONTINUE
      ELSE
*
*        Lower triangular storage: examine D from top to bottom.
*
         DO 20 INFO = 1, N
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   20    CONTINUE
      END IF
      INFO = 0
*
      IF( UPPER ) THEN
*
*        Compute inv(A) from the factorization A = U*D*U**H.
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         K = 1
   30    CONTINUE
*
*        If K > N, exit from loop.
*
         IF( K.GT.N )
     $      GO TO 50
*
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Invert the diagonal block.
*
            A( K, K ) = ONE / REAL( A( K, K ) )
*
*           Compute column K of the inverse.
*
            IF( K.GT.1 ) THEN
               CALL CCOPY( K-1, A( 1, K ), 1, WORK, 1 )
               CALL CHEMV( UPLO, K-1, -CONE, A, LDA, WORK, 1, ZERO,
     $                     A( 1, K ), 1 )
               A( K, K ) = A( K, K ) - REAL( CDOTC( K-1, WORK, 1, A( 1,
     $                     K ), 1 ) )
            END IF
            KSTEP = 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Invert the diagonal block.
*
            T = ABS( A( K, K+1 ) )
            AK = REAL( A( K, K ) ) / T
            AKP1 = REAL( A( K+1, K+1 ) ) / T
            AKKP1 = A( K, K+1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K, K ) = AKP1 / D
            A( K+1, K+1 ) = AK / D
            A( K, K+1 ) = -AKKP1 / D
*
*           Compute columns K and K+1 of the inverse.
*
            IF( K.GT.1 ) THEN
               CALL CCOPY( K-1, A( 1, K ), 1, WORK, 1 )
               CALL CHEMV( UPLO, K-1, -CONE, A, LDA, WORK, 1, ZERO,
     $                     A( 1, K ), 1 )
               A( K, K ) = A( K, K ) - REAL( CDOTC( K-1, WORK, 1, A( 1,
     $                     K ), 1 ) )
               A( K, K+1 ) = A( K, K+1 ) -
     $                       CDOTC( K-1, A( 1, K ), 1, A( 1, K+1 ), 1 )
               CALL CCOPY( K-1, A( 1, K+1 ), 1, WORK, 1 )
               CALL CHEMV( UPLO, K-1, -CONE, A, LDA, WORK, 1, ZERO,
     $                     A( 1, K+1 ), 1 )
               A( K+1, K+1 ) = A( K+1, K+1 ) -
     $                         REAL( CDOTC( K-1, WORK, 1, A( 1, K+1 ),
     $                         1 ) )
            END IF
            KSTEP = 2
         END IF
*
         KP = ABS( IPIV( K ) )
         IF( KP.NE.K ) THEN
*
*           Interchange rows and columns K and KP in the leading
*           submatrix A(1:k+1,1:k+1)
*
            CALL CSWAP( KP-1, A( 1, K ), 1, A( 1, KP ), 1 )
            DO 40 J = KP + 1, K - 1
               TEMP = CONJG( A( J, K ) )
               A( J, K ) = CONJG( A( KP, J ) )
               A( KP, J ) = TEMP
   40       CONTINUE
            A( KP, K ) = CONJG( A( KP, K ) )
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
            IF( KSTEP.EQ.2 ) THEN
               TEMP = A( K, K+1 )
               A( K, K+1 ) = A( KP, K+1 )
               A( KP, K+1 ) = TEMP
            END IF
         END IF
*
         K = K + KSTEP
         GO TO 30
   50    CONTINUE
*
      ELSE
*
*        Compute inv(A) from the factorization A = L*D*L**H.
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         K = N
   60    CONTINUE
*
*        If K < 1, exit from loop.
*
         IF( K.LT.1 )
     $      GO TO 80
*
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Invert the diagonal block.
*
            A( K, K ) = ONE / REAL( A( K, K ) )
*
*           Compute column K of the inverse.
*
            IF( K.LT.N ) THEN
               CALL CCOPY( N-K, A( K+1, K ), 1, WORK, 1 )
               CALL CHEMV( UPLO, N-K, -CONE, A( K+1, K+1 ), LDA, WORK,
     $                     1, ZERO, A( K+1, K ), 1 )
               A( K, K ) = A( K, K ) - REAL( CDOTC( N-K, WORK, 1,
     $                     A( K+1, K ), 1 ) )
            END IF
            KSTEP = 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Invert the diagonal block.
*
            T = ABS( A( K, K-1 ) )
            AK = REAL( A( K-1, K-1 ) ) / T
            AKP1 = REAL( A( K, K ) ) / T
            AKKP1 = A( K, K-1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K-1, K-1 ) = AKP1 / D
            A( K, K ) = AK / D
            A( K, K-1 ) = -AKKP1 / D
*
*           Compute columns K-1 and K of the inverse.
*
            IF( K.LT.N ) THEN
               CALL CCOPY( N-K, A( K+1, K ), 1, WORK, 1 )
               CALL CHEMV( UPLO, N-K, -CONE, A( K+1, K+1 ), LDA, WORK,
     $                     1, ZERO, A( K+1, K ), 1 )
               A( K, K ) = A( K, K ) - REAL( CDOTC( N-K, WORK, 1,
     $                     A( K+1, K ), 1 ) )
               A( K, K-1 ) = A( K, K-1 ) -
     $                       CDOTC( N-K, A( K+1, K ), 1, A( K+1, K-1 ),
     $                       1 )
               CALL CCOPY( N-K, A( K+1, K-1 ), 1, WORK, 1 )
               CALL CHEMV( UPLO, N-K, -CONE, A( K+1, K+1 ), LDA, WORK,
     $                     1, ZERO, A( K+1, K-1 ), 1 )
               A( K-1, K-1 ) = A( K-1, K-1 ) -
     $                         REAL( CDOTC( N-K, WORK, 1, A( K+1, K-1 ),
     $                         1 ) )
            END IF
            KSTEP = 2
         END IF
*
         KP = ABS( IPIV( K ) )
         IF( KP.NE.K ) THEN
*
*           Interchange rows and columns K and KP in the trailing
*           submatrix A(k-1:n,k-1:n)
*
            IF( KP.LT.N )
     $         CALL CSWAP( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 )
            DO 70 J = K + 1, KP - 1
               TEMP = CONJG( A( J, K ) )
               A( J, K ) = CONJG( A( KP, J ) )
               A( KP, J ) = TEMP
   70       CONTINUE
            A( KP, K ) = CONJG( A( KP, K ) )
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
            IF( KSTEP.EQ.2 ) THEN
               TEMP = A( K, K-1 )
               A( K, K-1 ) = A( KP, K-1 )
               A( KP, K-1 ) = TEMP
            END IF
         END IF
*
         K = K - KSTEP
         GO TO 60
   80    CONTINUE
      END IF
*
      RETURN
*
*     End of CHETRI
*
      END
