# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/clarzb.f"
      SUBROUTINE CLARZB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, V,
     $                   LDV, T, LDT, C, LDC, WORK, LDWORK )
*
*  -- LAPACK routine (version 3.3.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- April 2011                                                      --
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, L, LDC, LDT, LDV, LDWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX            C( LDC, * ), T( LDT, * ), V( LDV, * ),
     $                   WORK( LDWORK, * )
*     ..
*
*  Purpose
*  =======
*
*  CLARZB applies a complex block reflector H or its transpose H**H
*  to a complex distributed M-by-N  C from the left or the right.
*
*  Currently, only STOREV = 'R' and DIRECT = 'B' are supported.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply H or H**H from the Left
*          = 'R': apply H or H**H from the Right
*
*  TRANS   (input) CHARACTER*1
*          = 'N': apply H (No transpose)
*          = 'C': apply H**H (Conjugate transpose)
*
*  DIRECT  (input) CHARACTER*1
*          Indicates how H is formed from a product of elementary
*          reflectors
*          = 'F': H = H(1) H(2) . . . H(k) (Forward, not supported yet)
*          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (input) CHARACTER*1
*          Indicates how the vectors which define the elementary
*          reflectors are stored:
*          = 'C': Columnwise                        (not supported yet)
*          = 'R': Rowwise
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  K       (input) INTEGER
*          The order of the matrix T (= the number of elementary
*          reflectors whose product defines the block reflector).
*
*  L       (input) INTEGER
*          The number of columns of the matrix V containing the
*          meaningful part of the Householder reflectors.
*          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.
*
*  V       (input) COMPLEX array, dimension (LDV,NV).
*          If STOREV = 'C', NV = K; if STOREV = 'R', NV = L.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V.
*          If STOREV = 'C', LDV >= L; if STOREV = 'R', LDV >= K.
*
*  T       (input) COMPLEX array, dimension (LDT,K)
*          The triangular K-by-K matrix T in the representation of the
*          block reflector.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= K.
*
*  C       (input/output) COMPLEX array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by H*C or H**H*C or C*H or C*H**H.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) COMPLEX array, dimension (LDWORK,K)
*
*  LDWORK  (input) INTEGER
*          The leading dimension of the array WORK.
*          If SIDE = 'L', LDWORK >= max(1,N);
*          if SIDE = 'R', LDWORK >= max(1,M).
*
*  Further Details
*  ===============
*
*  Based on contributions by
*    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, INFO, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CCOPY, CGEMM, CLACGV, CTRMM, XERBLA
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
*     Check for currently supported options
*
      INFO = 0
      IF( .NOT.LSAME( DIRECT, 'B' ) ) THEN
         INFO = -3
      ELSE IF( .NOT.LSAME( STOREV, 'R' ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CLARZB', -INFO )
         RETURN
      END IF
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'C'
      ELSE
         TRANST = 'N'
      END IF
*
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  H * C  or  H**H * C
*
*        W( 1:n, 1:k ) = C( 1:k, 1:n )**H
*
         DO 10 J = 1, K
            CALL CCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
   10    CONTINUE
*
*        W( 1:n, 1:k ) = W( 1:n, 1:k ) + ...
*                        C( m-l+1:m, 1:n )**H * V( 1:k, 1:l )**T
*
         IF( L.GT.0 )
     $      CALL CGEMM( 'Transpose', 'Conjugate transpose', N, K, L,
     $                  ONE, C( M-L+1, 1 ), LDC, V, LDV, ONE, WORK,
     $                  LDWORK )
*
*        W( 1:n, 1:k ) = W( 1:n, 1:k ) * T**T  or  W( 1:m, 1:k ) * T
*
         CALL CTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, ONE, T,
     $               LDT, WORK, LDWORK )
*
*        C( 1:k, 1:n ) = C( 1:k, 1:n ) - W( 1:n, 1:k )**H
*
         DO 30 J = 1, N
            DO 20 I = 1, K
               C( I, J ) = C( I, J ) - WORK( J, I )
   20       CONTINUE
   30    CONTINUE
*
*        C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
*                            V( 1:k, 1:l )**H * W( 1:n, 1:k )**H
*
         IF( L.GT.0 )
     $      CALL CGEMM( 'Transpose', 'Transpose', L, N, K, -ONE, V, LDV,
     $                  WORK, LDWORK, ONE, C( M-L+1, 1 ), LDC )
*
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*        Form  C * H  or  C * H**H
*
*        W( 1:m, 1:k ) = C( 1:m, 1:k )
*
         DO 40 J = 1, K
            CALL CCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40    CONTINUE
*
*        W( 1:m, 1:k ) = W( 1:m, 1:k ) + ...
*                        C( 1:m, n-l+1:n ) * V( 1:k, 1:l )**H
*
         IF( L.GT.0 )
     $      CALL CGEMM( 'No transpose', 'Transpose', M, K, L, ONE,
     $                  C( 1, N-L+1 ), LDC, V, LDV, ONE, WORK, LDWORK )
*
*        W( 1:m, 1:k ) = W( 1:m, 1:k ) * conjg( T )  or
*                        W( 1:m, 1:k ) * T**H
*
         DO 50 J = 1, K
            CALL CLACGV( K-J+1, T( J, J ), 1 )
   50    CONTINUE
         CALL CTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, ONE, T,
     $               LDT, WORK, LDWORK )
         DO 60 J = 1, K
            CALL CLACGV( K-J+1, T( J, J ), 1 )
   60    CONTINUE
*
*        C( 1:m, 1:k ) = C( 1:m, 1:k ) - W( 1:m, 1:k )
*
         DO 80 J = 1, K
            DO 70 I = 1, M
               C( I, J ) = C( I, J ) - WORK( I, J )
   70       CONTINUE
   80    CONTINUE
*
*        C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
*                            W( 1:m, 1:k ) * conjg( V( 1:k, 1:l ) )
*
         DO 90 J = 1, L
            CALL CLACGV( K, V( 1, J ), 1 )
   90    CONTINUE
         IF( L.GT.0 )
     $      CALL CGEMM( 'No transpose', 'No transpose', M, L, K, -ONE,
     $                  WORK, LDWORK, V, LDV, ONE, C( 1, N-L+1 ), LDC )
         DO 100 J = 1, L
            CALL CLACGV( K, V( 1, J ), 1 )
  100    CONTINUE
*
      END IF
*
      RETURN
*
*     End of CLARZB
*
      END
