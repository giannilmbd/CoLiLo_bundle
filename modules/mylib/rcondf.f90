!     Last change:  GL   14 Nov 2006    5:29 pm
subroutine rcondf(A_,mm,rc)


implicit none

DOUBLE PRECISION, DIMENSION(mm,mm), INTENT(IN) ::A_
DOUBLE PRECISION, DIMENSION(mm,mm)::A
INTEGER*8, INTENT(IN) ::mm
double precision, INTENT(OUT) :: rc
DOUBLE PRECISION, DIMENSION(mm*4) :: work
INTEGER*8, DIMENSION(mm):: iwork
DOUBLE PRECISION, EXTERNAL :: Dlange
DOUBLE precision :: na
INTEGER*8 :: info ,ipiv(mm)
A=A_;
!first compute the norm of A
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER!1
!          Specifies the value to be returned in DLANGE as described
!          above.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.  When M = 0,
!          DLANGE is set to zero.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.  When N = 0,
!          DLANGE is set to zero.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The m by n matrix A.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(M,1).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
!          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
!          referenced.
!
! =====================================================================

na=DLANGE('1', mm, mm, A, mm, WORK )

!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P!L!U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================

call DGETRF( mm, mm, A, mm, IPIV, INFO )

!  Arguments
!  =========
!
!  NORM    (input) CHARACTER!1
!          Specifies whether the 1-norm condition number or the
!          infinity-norm condition number is required:
!          = '1' or 'O':  1-norm;
!          = 'I':         Infinity-norm.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The factors L and U from the factorization A = P!L!U
!          as computed by DGETRF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  ANORM   (input) DOUBLE PRECISION
!          If NORM = '1' or 'O', the 1-norm of the original matrix A.
!          If NORM = 'I', the infinity-norm of the original matrix A.
!
!  RCOND   (output) DOUBLE PRECISION
!          The reciprocal of the condition number of the matrix A,
!          computed as RCOND = 1/(norm(A) ! norm(inv(A))).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (4!N)
!
!  IWORK   (workspace) INTEGER array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================


call DGECON( '1', mm, A, mm, na,rc, WORK, IWORK,INFO )

end subroutine
