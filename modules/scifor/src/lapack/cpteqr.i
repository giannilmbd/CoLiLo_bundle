# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/cpteqr.f"
      SUBROUTINE CPTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
*     ..
*     .. Array Arguments ..
      REAL               D( * ), E( * ), WORK( * )
      COMPLEX            Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  CPTEQR computes all eigenvalues and, optionally, eigenvectors of a
*  symmetric positive definite tridiagonal matrix by first factoring the
*  matrix using SPTTRF and then calling CBDSQR to compute the singular
*  values of the bidiagonal factor.
*
*  This routine computes the eigenvalues of the positive definite
*  tridiagonal matrix to high relative accuracy.  This means that if the
*  eigenvalues range over many orders of magnitude in size, then the
*  small eigenvalues and corresponding eigenvectors will be computed
*  more accurately than, for example, with the standard QR method.
*
*  The eigenvectors of a full or band positive definite Hermitian matrix
*  can also be found if CHETRD, CHPTRD, or CHBTRD has been used to
*  reduce this matrix to tridiagonal form.  (The reduction to
*  tridiagonal form, however, may preclude the possibility of obtaining
*  high relative accuracy in the small eigenvalues of the original
*  matrix, if these eigenvalues range over many orders of magnitude.)
*
*  Arguments
*  =========
*
*  COMPZ   (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only.
*          = 'V':  Compute eigenvectors of original Hermitian
*                  matrix also.  Array Z contains the unitary matrix
*                  used to reduce the original matrix to tridiagonal
*                  form.
*          = 'I':  Compute eigenvectors of tridiagonal matrix also.
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) REAL array, dimension (N)
*          On entry, the n diagonal elements of the tridiagonal matrix.
*          On normal exit, D contains the eigenvalues, in descending
*          order.
*
*  E       (input/output) REAL array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix.
*          On exit, E has been destroyed.
*
*  Z       (input/output) COMPLEX array, dimension (LDZ, N)
*          On entry, if COMPZ = 'V', the unitary matrix used in the
*          reduction to tridiagonal form.
*          On exit, if COMPZ = 'V', the orthonormal eigenvectors of the
*          original Hermitian matrix;
*          if COMPZ = 'I', the orthonormal eigenvectors of the
*          tridiagonal matrix.
*          If INFO > 0 on exit, Z contains the eigenvectors associated
*          with only the stored eigenvalues.
*          If  COMPZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          COMPZ = 'V' or 'I', LDZ >= max(1,N).
*
*  WORK    (workspace) REAL array, dimension (4*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = i, and i is:
*                <= N  the Cholesky factorization of the matrix could
*                      not be performed because the i-th principal minor
*                      was not positive definite.
*                > N   the SVD algorithm failed to converge;
*                      if INFO = N+i, i off-diagonal elements of the
*                      bidiagonal factor did not converge to zero.
*
*  ====================================================================
*
*     .. Parameters ..
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ),
     $                   CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CBDSQR, CLASET, SPTTRF, XERBLA
*     ..
*     .. Local Arrays ..
      COMPLEX            C( 1, 1 ), VT( 1, 1 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICOMPZ, NRU
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
     $         N ) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CPTEQR', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.GT.0 )
     $      Z( 1, 1 ) = CONE
         RETURN
      END IF
      IF( ICOMPZ.EQ.2 )
     $   CALL CLASET( 'Full', N, N, CZERO, CONE, Z, LDZ )
*
*     Call SPTTRF to factor the matrix.
*
      CALL SPTTRF( N, D, E, INFO )
      IF( INFO.NE.0 )
     $   RETURN
      DO 10 I = 1, N
         D( I ) = SQRT( D( I ) )
   10 CONTINUE
      DO 20 I = 1, N - 1
         E( I ) = E( I )*D( I )
   20 CONTINUE
*
*     Call CBDSQR to compute the singular values/vectors of the
*     bidiagonal factor.
*
      IF( ICOMPZ.GT.0 ) THEN
         NRU = N
      ELSE
         NRU = 0
      END IF
      CALL CBDSQR( 'Lower', N, 0, NRU, 0, D, E, VT, 1, Z, LDZ, C, 1,
     $             WORK, INFO )
*
*     Square the singular values.
*
      IF( INFO.EQ.0 ) THEN
         DO 30 I = 1, N
            D( I ) = D( I )*D( I )
   30    CONTINUE
      ELSE
         INFO = N + INFO
      END IF
*
      RETURN
*
*     End of CPTEQR
*
      END
