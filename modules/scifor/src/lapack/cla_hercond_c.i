# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/cla_hercond_c.f"
      REAL FUNCTION CLA_HERCOND_C( UPLO, N, A, LDA, AF, LDAF, IPIV, C,
     $                             CAPPLY, INFO, WORK, RWORK )
*
*     -- LAPACK routine (version 3.2.1)                                 --
*     -- Contributed by James Demmel, Deaglan Halligan, Yozo Hida and --
*     -- Jason Riedy of Univ. of California Berkeley.                 --
*     -- April 2009                                                   --
*
*     -- LAPACK is a software package provided by Univ. of Tennessee, --
*     -- Univ. of California Berkeley and NAG Ltd.                    --
*
      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      LOGICAL            CAPPLY
      INTEGER            N, LDA, LDAF, INFO
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * ), AF( LDAF, * ), WORK( * )
      REAL               C ( * ), RWORK( * )
*     ..
*
*  Purpose
*  =======
*
*     CLA_HERCOND_C computes the infinity norm condition number of
*     op(A) * inv(diag(C)) where C is a REAL vector.
*
*  Arguments
*  =========
*
*     UPLO    (input) CHARACTER*1
*       = 'U':  Upper triangle of A is stored;
*       = 'L':  Lower triangle of A is stored.
*
*     N       (input) INTEGER
*     The number of linear equations, i.e., the order of the
*     matrix A.  N >= 0.
*
*     A       (input) COMPLEX array, dimension (LDA,N)
*     On entry, the N-by-N matrix A
*
*     LDA     (input) INTEGER
*     The leading dimension of the array A.  LDA >= max(1,N).
*
*     AF      (input) COMPLEX array, dimension (LDAF,N)
*     The block diagonal matrix D and the multipliers used to
*     obtain the factor U or L as computed by CHETRF.
*
*     LDAF    (input) INTEGER
*     The leading dimension of the array AF.  LDAF >= max(1,N).
*
*     IPIV    (input) INTEGER array, dimension (N)
*     Details of the interchanges and the block structure of D
*     as determined by CHETRF.
*
*     C       (input) REAL array, dimension (N)
*     The vector C in the formula op(A) * inv(diag(C)).
*
*     CAPPLY  (input) LOGICAL
*     If .TRUE. then access the vector C in the formula above.
*
*     INFO    (output) INTEGER
*       = 0:  Successful exit.
*     i > 0:  The ith argument is invalid.
*
*     WORK    (input) COMPLEX array, dimension (2*N).
*     Workspace.
*
*     RWORK   (input) REAL array, dimension (N).
*     Workspace.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            KASE, I, J
      REAL               AINVNM, ANORM, TMP
      LOGICAL            UP
      COMPLEX            ZDUM
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLACN2, CHETRS, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
      CLA_HERCOND_C = 0.0E+0
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CLA_HERCOND_C', -INFO )
         RETURN
      END IF
      UP = .FALSE.
      IF ( LSAME( UPLO, 'U' ) ) UP = .TRUE.
*
*     Compute norm of op(A)*op2(C).
*
      ANORM = 0.0E+0
      IF ( UP ) THEN
         DO I = 1, N
            TMP = 0.0E+0
            IF ( CAPPLY ) THEN
               DO J = 1, I
                  TMP = TMP + CABS1( A( J, I ) ) / C( J )
               END DO
               DO J = I+1, N
                  TMP = TMP + CABS1( A( I, J ) ) / C( J )
               END DO
            ELSE
               DO J = 1, I
                  TMP = TMP + CABS1( A( J, I ) )
               END DO
               DO J = I+1, N
                  TMP = TMP + CABS1( A( I, J ) )
               END DO
            END IF
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      ELSE
         DO I = 1, N
            TMP = 0.0E+0
            IF ( CAPPLY ) THEN
               DO J = 1, I
                  TMP = TMP + CABS1( A( I, J ) ) / C( J )
               END DO
               DO J = I+1, N
                  TMP = TMP + CABS1( A( J, I ) ) / C( J )
               END DO
            ELSE
               DO J = 1, I
                  TMP = TMP + CABS1( A( I, J ) )
               END DO
               DO J = I+1, N
                  TMP = TMP + CABS1( A( J, I ) )
               END DO
            END IF
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         END DO
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 ) THEN
         CLA_HERCOND_C = 1.0E+0
         RETURN
      ELSE IF( ANORM .EQ. 0.0E+0 ) THEN
         RETURN
      END IF
*
*     Estimate the norm of inv(op(A)).
*
      AINVNM = 0.0E+0
*
      KASE = 0
   10 CONTINUE
      CALL CLACN2( N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.2 ) THEN
*
*           Multiply by R.
*
            DO I = 1, N
               WORK( I ) = WORK( I ) * RWORK( I )
            END DO
*
            IF ( UP ) THEN
               CALL CHETRS( 'U', N, 1, AF, LDAF, IPIV,
     $            WORK, N, INFO )
            ELSE
               CALL CHETRS( 'L', N, 1, AF, LDAF, IPIV,
     $            WORK, N, INFO )
            ENDIF
*
*           Multiply by inv(C).
*
            IF ( CAPPLY ) THEN
               DO I = 1, N
                  WORK( I ) = WORK( I ) * C( I )
               END DO
            END IF
         ELSE
*
*           Multiply by inv(C**H).
*
            IF ( CAPPLY ) THEN
               DO I = 1, N
                  WORK( I ) = WORK( I ) * C( I )
               END DO
            END IF
*
            IF ( UP ) THEN
               CALL CHETRS( 'U', N, 1, AF, LDAF, IPIV,
     $            WORK, N, INFO )
            ELSE
               CALL CHETRS( 'L', N, 1, AF, LDAF, IPIV,
     $            WORK, N, INFO )
            END IF
*
*           Multiply by R.
*
            DO I = 1, N
               WORK( I ) = WORK( I ) * RWORK( I )
            END DO
         END IF
         GO TO 10
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM .NE. 0.0E+0 )
     $   CLA_HERCOND_C = 1.0E+0 / AINVNM
*
      RETURN
*
      END
