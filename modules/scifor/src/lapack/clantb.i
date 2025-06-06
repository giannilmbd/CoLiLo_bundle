# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/clantb.f"
      REAL             FUNCTION CLANTB( NORM, UPLO, DIAG, N, K, AB,
     $                 LDAB, WORK )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, NORM, UPLO
      INTEGER            K, LDAB, N
*     ..
*     .. Array Arguments ..
      REAL               WORK( * )
      COMPLEX            AB( LDAB, * )
*     ..
*
*  Purpose
*  =======
*
*  CLANTB  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the element of  largest absolute value  of an
*  n by n triangular band matrix A,  with ( k + 1 ) diagonals.
*
*  Description
*  ===========
*
*  CLANTB returns the value
*
*     CLANTB = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in CLANTB as described
*          above.
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the matrix A is upper or lower triangular.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  DIAG    (input) CHARACTER*1
*          Specifies whether or not the matrix A is unit triangular.
*          = 'N':  Non-unit triangular
*          = 'U':  Unit triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.  When N = 0, CLANTB is
*          set to zero.
*
*  K       (input) INTEGER
*          The number of super-diagonals of the matrix A if UPLO = 'U',
*          or the number of sub-diagonals of the matrix A if UPLO = 'L'.
*          K >= 0.
*
*  AB      (input) COMPLEX array, dimension (LDAB,N)
*          The upper or lower triangular band matrix A, stored in the
*          first k+1 rows of AB.  The j-th column of A is stored
*          in the j-th column of the array AB as follows:
*          if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k).
*          Note that when DIAG = 'U', the elements of the array AB
*          corresponding to the diagonal elements of the matrix A are
*          not referenced, but are assumed to be one.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= K+1.
*
*  WORK    (workspace) REAL array, dimension (MAX(1,LWORK)),
*          where LWORK >= N when NORM = 'I'; otherwise, WORK is not
*          referenced.
*
* =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UDIAG
      INTEGER            I, J, L
      REAL               SCALE, SUM, VALUE
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLASSQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         IF( LSAME( DIAG, 'U' ) ) THEN
            VALUE = ONE
            IF( LSAME( UPLO, 'U' ) ) THEN
               DO 20 J = 1, N
                  DO 10 I = MAX( K+2-J, 1 ), K
                     VALUE = MAX( VALUE, ABS( AB( I, J ) ) )
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40 J = 1, N
                  DO 30 I = 2, MIN( N+1-J, K+1 )
                     VALUE = MAX( VALUE, ABS( AB( I, J ) ) )
   30             CONTINUE
   40          CONTINUE
            END IF
         ELSE
            VALUE = ZERO
            IF( LSAME( UPLO, 'U' ) ) THEN
               DO 60 J = 1, N
                  DO 50 I = MAX( K+2-J, 1 ), K + 1
                     VALUE = MAX( VALUE, ABS( AB( I, J ) ) )
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80 J = 1, N
                  DO 70 I = 1, MIN( N+1-J, K+1 )
                     VALUE = MAX( VALUE, ABS( AB( I, J ) ) )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
*
*        Find norm1(A).
*
         VALUE = ZERO
         UDIAG = LSAME( DIAG, 'U' )
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 1, N
               IF( UDIAG ) THEN
                  SUM = ONE
                  DO 90 I = MAX( K+2-J, 1 ), K
                     SUM = SUM + ABS( AB( I, J ) )
   90             CONTINUE
               ELSE
                  SUM = ZERO
                  DO 100 I = MAX( K+2-J, 1 ), K + 1
                     SUM = SUM + ABS( AB( I, J ) )
  100             CONTINUE
               END IF
               VALUE = MAX( VALUE, SUM )
  110       CONTINUE
         ELSE
            DO 140 J = 1, N
               IF( UDIAG ) THEN
                  SUM = ONE
                  DO 120 I = 2, MIN( N+1-J, K+1 )
                     SUM = SUM + ABS( AB( I, J ) )
  120             CONTINUE
               ELSE
                  SUM = ZERO
                  DO 130 I = 1, MIN( N+1-J, K+1 )
                     SUM = SUM + ABS( AB( I, J ) )
  130             CONTINUE
               END IF
               VALUE = MAX( VALUE, SUM )
  140       CONTINUE
         END IF
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
*
*        Find normI(A).
*
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            IF( LSAME( DIAG, 'U' ) ) THEN
               DO 150 I = 1, N
                  WORK( I ) = ONE
  150          CONTINUE
               DO 170 J = 1, N
                  L = K + 1 - J
                  DO 160 I = MAX( 1, J-K ), J - 1
                     WORK( I ) = WORK( I ) + ABS( AB( L+I, J ) )
  160             CONTINUE
  170          CONTINUE
            ELSE
               DO 180 I = 1, N
                  WORK( I ) = ZERO
  180          CONTINUE
               DO 200 J = 1, N
                  L = K + 1 - J
                  DO 190 I = MAX( 1, J-K ), J
                     WORK( I ) = WORK( I ) + ABS( AB( L+I, J ) )
  190             CONTINUE
  200          CONTINUE
            END IF
         ELSE
            IF( LSAME( DIAG, 'U' ) ) THEN
               DO 210 I = 1, N
                  WORK( I ) = ONE
  210          CONTINUE
               DO 230 J = 1, N
                  L = 1 - J
                  DO 220 I = J + 1, MIN( N, J+K )
                     WORK( I ) = WORK( I ) + ABS( AB( L+I, J ) )
  220             CONTINUE
  230          CONTINUE
            ELSE
               DO 240 I = 1, N
                  WORK( I ) = ZERO
  240          CONTINUE
               DO 260 J = 1, N
                  L = 1 - J
                  DO 250 I = J, MIN( N, J+K )
                     WORK( I ) = WORK( I ) + ABS( AB( L+I, J ) )
  250             CONTINUE
  260          CONTINUE
            END IF
         END IF
         DO 270 I = 1, N
            VALUE = MAX( VALUE, WORK( I ) )
  270    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         IF( LSAME( UPLO, 'U' ) ) THEN
            IF( LSAME( DIAG, 'U' ) ) THEN
               SCALE = ONE
               SUM = N
               IF( K.GT.0 ) THEN
                  DO 280 J = 2, N
                     CALL CLASSQ( MIN( J-1, K ),
     $                            AB( MAX( K+2-J, 1 ), J ), 1, SCALE,
     $                            SUM )
  280             CONTINUE
               END IF
            ELSE
               SCALE = ZERO
               SUM = ONE
               DO 290 J = 1, N
                  CALL CLASSQ( MIN( J, K+1 ), AB( MAX( K+2-J, 1 ), J ),
     $                         1, SCALE, SUM )
  290          CONTINUE
            END IF
         ELSE
            IF( LSAME( DIAG, 'U' ) ) THEN
               SCALE = ONE
               SUM = N
               IF( K.GT.0 ) THEN
                  DO 300 J = 1, N - 1
                     CALL CLASSQ( MIN( N-J, K ), AB( 2, J ), 1, SCALE,
     $                            SUM )
  300             CONTINUE
               END IF
            ELSE
               SCALE = ZERO
               SUM = ONE
               DO 310 J = 1, N
                  CALL CLASSQ( MIN( N-J+1, K+1 ), AB( 1, J ), 1, SCALE,
     $                         SUM )
  310          CONTINUE
            END IF
         END IF
         VALUE = SCALE*SQRT( SUM )
      END IF
*
      CLANTB = VALUE
      RETURN
*
*     End of CLANTB
*
      END
