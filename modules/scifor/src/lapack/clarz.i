# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/clarz.f"
      SUBROUTINE CLARZ( SIDE, M, N, L, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK routine (version 3.3.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- April 2011                                                      --
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, L, LDC, M, N
      COMPLEX            TAU
*     ..
*     .. Array Arguments ..
      COMPLEX            C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CLARZ applies a complex elementary reflector H to a complex
*  M-by-N matrix C, from either the left or the right. H is represented
*  in the form
*
*        H = I - tau * v * v**H
*
*  where tau is a complex scalar and v is a complex vector.
*
*  If tau = 0, then H is taken to be the unit matrix.
*
*  To apply H**H (the conjugate transpose of H), supply conjg(tau) instead
*  tau.
*
*  H is a product of k elementary reflectors as returned by CTZRZF.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': form  H * C
*          = 'R': form  C * H
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  L       (input) INTEGER
*          The number of entries of the vector V containing
*          the meaningful part of the Householder vectors.
*          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.
*
*  V       (input) COMPLEX array, dimension (1+(L-1)*abs(INCV))
*          The vector v in the representation of H as returned by
*          CTZRZF. V is not used if TAU = 0.
*
*  INCV    (input) INTEGER
*          The increment between elements of v. INCV <> 0.
*
*  TAU     (input) COMPLEX
*          The value tau in the representation of H.
*
*  C       (input/output) COMPLEX array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
*          or C * H if SIDE = 'R'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) COMPLEX array, dimension
*                         (N) if SIDE = 'L'
*                      or (M) if SIDE = 'R'
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
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),
     $                   ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. External Subroutines ..
      EXTERNAL           CAXPY, CCOPY, CGEMV, CGERC, CGERU, CLACGV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  H * C
*
         IF( TAU.NE.ZERO ) THEN
*
*           w( 1:n ) = conjg( C( 1, 1:n ) )
*
            CALL CCOPY( N, C, LDC, WORK, 1 )
            CALL CLACGV( N, WORK, 1 )
*
*           w( 1:n ) = conjg( w( 1:n ) + C( m-l+1:m, 1:n )**H * v( 1:l ) )
*
            CALL CGEMV( 'Conjugate transpose', L, N, ONE, C( M-L+1, 1 ),
     $                  LDC, V, INCV, ONE, WORK, 1 )
            CALL CLACGV( N, WORK, 1 )
*
*           C( 1, 1:n ) = C( 1, 1:n ) - tau * w( 1:n )
*
            CALL CAXPY( N, -TAU, WORK, 1, C, LDC )
*
*           C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
*                               tau * v( 1:l ) * w( 1:n )**H
*
            CALL CGERU( L, N, -TAU, V, INCV, WORK, 1, C( M-L+1, 1 ),
     $                  LDC )
         END IF
*
      ELSE
*
*        Form  C * H
*
         IF( TAU.NE.ZERO ) THEN
*
*           w( 1:m ) = C( 1:m, 1 )
*
            CALL CCOPY( M, C, 1, WORK, 1 )
*
*           w( 1:m ) = w( 1:m ) + C( 1:m, n-l+1:n, 1:n ) * v( 1:l )
*
            CALL CGEMV( 'No transpose', M, L, ONE, C( 1, N-L+1 ), LDC,
     $                  V, INCV, ONE, WORK, 1 )
*
*           C( 1:m, 1 ) = C( 1:m, 1 ) - tau * w( 1:m )
*
            CALL CAXPY( M, -TAU, WORK, 1, C, 1 )
*
*           C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
*                               tau * w( 1:m ) * v( 1:l )**H
*
            CALL CGERC( M, L, -TAU, WORK, 1, V, INCV, C( 1, N-L+1 ),
     $                  LDC )
*
         END IF
*
      END IF
*
      RETURN
*
*     End of CLARZ
*
      END
