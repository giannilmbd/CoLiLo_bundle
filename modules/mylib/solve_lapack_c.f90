
 SUBROUTINE solve_lapack_c(Ain,LDAin,nin,Bin,nbin,Xout,info)

implicit none

! input/output of main subroutine
integer*8 :: LDAin,nin,nbin
complex*16 :: Ain(LDAin,nin),Bin(LDAin,nbin)

complex*16:: Xout(nin,nbin)
! inputs of the factorization
! character :: UPLO
integer*8 :: N,LDA,INFO
integer*8, dimension(Nin) :: IPIV
complex*16 :: A(LDAin,Nin),B(LDAin,nbin)
! double precision, allocatable, dimension(:) :: WORK
integer*8 :: error
! external zsytrf,dsytrs

N=nin
LDA=LDAin
A=Ain;
B=Bin;

Xout=0.0;

INFO=0;


! factorization: use this if use ZGETRS
!       call ZGETRF( LDA, N, A, LDA, IPIV, INFO )


! I used to use this before using Matlab
! call ZGETRS( 'N', N, nbin, A, LDA, IPIV, B, LDA, INFO )





! the following is used by Matlab

! PURPOSE
!      ZGESV computes the solution to a complex system of linear
!      equations
!         A * X = B, where A is an N-by-N matrix and X and B are
!      N-by-NRHS matrices.

!      The LU decomposition with partial pivoting and row inter-
!      changes is used to factor A as
!         A = P * L * U,
!      where P is a permutation matrix, L is unit lower triangular,
!      and U is upper triangular.  The factored form of A is then
!      used to solve the system of equations A * X = B.
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of linear equations, i.e., the order of the
!          matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the N-by-N coefficient matrix A.
!          On exit, the factors L and U from the factorization
!          A = P!L!U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (output) INTEGER array, dimension (N)
!          The pivot indices that define the permutation matrix P;
!          row i of the matrix was interchanged with row IPIV(i).
!
!  B       (input/output) COMPLEX!16 array, dimension (LDB,NRHS)
!          On entry, the N-by-NRHS matrix of right hand side matrix B.
!          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
!                has been completed, but the factor U is exactly
!                singular, so the solution could not be computed.
!
!  =====================================================================
!  ZGESV - compute the solution to a complex system of linear
!      equations  A * X = B,
call ZGESV( N, nbin, A, LDA, IPIV, B, N, INFO )

Xout=B;

end subroutine
