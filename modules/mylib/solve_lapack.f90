
 SUBROUTINE solve_lapack(Ain,LDAin,nin,Bin,nbin,Xout)

implicit none

! input/output of main subroutine
integer, intent(in) :: LDAin,nin,nbin
double precision, intent(in) :: Ain(LDAin,nin),Bin(LDAin,nbin)

double precision, intent(out) :: Xout(nin,nbin)
! inputs of the factorization
character :: UPLO
integer :: N,LDA,INFO
integer, dimension(Nin) :: IPIV
double precision :: A(LDAin,Nin),B(LDAin,nbin)
! double precision, allocatable, dimension(:) :: WORK
integer :: error
! external dsytrf,dsytrs

N=nin
LDA=LDAin
A=Ain;
B=Bin;

Xout=0.0d0;
! UPLO='U';
INFO=-1;
! allocate(work(1),stat=error)
! work(1)=1;
! Lwork=-1;

! factorization
      call     DGETRF( LDA, N, A, LDA, IPIV, INFO )
!       call DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
! Lwork=aint(work(1));
! deallocate(work)
! allocate(work(Lwork),stat=error)

!  call DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )

! deallocate(work)
! solution
!      DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
call   DGETRS( 'N', N, nbin, A, LDA, IPIV, B, LDA, INFO )
! call DSYTRS( UPLO, N, nbin, A, LDA, IPIV, B, LDA, INFO )

Xout=B;
! write(*,*) 'info',info
end subroutine
