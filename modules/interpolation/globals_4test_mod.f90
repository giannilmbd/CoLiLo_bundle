
module globals_4test_mod
    implicit none
    integer,parameter :: ngrid_=10  
    integer,parameter :: ord=2
    integer, parameter:: num_coeffs=(ord*ord)
    integer,parameter :: whichinv=3
    logical,parameter :: cheb=.true.
    integer, parameter :: iters=2000,iters_p=2000, iters_k=2000
    type :: params
    integer :: num_coeffs=ord*ord
    real*8 :: new_b_p(0:num_coeffs-1)
    real*8 :: new_b_XX(0:num_coeffs-1,2)
    real*8 :: new_b_kappa(0:num_coeffs-1)
    real*8 :: OnGrid(0:ngrid_-1,0:num_coeffs-1)
    real*8 :: Oinv(0:num_coeffs-1,0:ngrid_-1)
    real*8 :: Oridge(0:num_coeffs-1,0:ngrid_-1)
    real*8 :: LS(0:num_coeffs-1,0:ngrid_-1)
    real*8 :: Intg(0:ngrid_-1,0:num_coeffs-1)
    real*8 :: D(0:ngrid_-1)
    real*8 :: Ds(0:ngrid_-1)
    end type
    real*8,dimension(0:ngrid_-1):: D,Ds,D_st,Ds_st,D_stst,Ds_stst
    real*8 :: B1(0:ngrid_-1,0:ord-1), B2(0:ngrid_-1,0:ord-1)!,B2d(1:ngrid_,1:ord),B1d(1:ngrid_,1:ord)
    real*8:: A(0:ngrid_-1,0:num_coeffs-1)
    real*8 :: SS(0:num_coeffs-1,0:num_coeffs-1),Ssvd(0:num_coeffs-1)
    real*8 :: Usvd(0:ngrid_-1,0:num_coeffs-1),Vsvd(0:num_coeffs-1,0:num_coeffs-1)
    real*8, allocatable:: Ur(:,:),Vr(:,:),SSr(:,:)
    real*8 :: Xout(0:num_coeffs-1,0:num_coeffs-1),eye(0:num_coeffs-1,0:num_coeffs-1)
    real*8 :: param_p0(0:num_coeffs-1),param_X0(0:num_coeffs-1)
    real*8, dimension(0:ngrid_-1,0:num_coeffs-1) :: Ones


    contains
    function policy_function(vect,matx) result(new_vect)
        real*8, intent(in) :: vect(0:num_coeffs-1) 
        real*8, intent(in) :: matx(0:ngrid_-1,0:num_coeffs-1)
        real*8 :: new_vect(0:ngrid_-1)

        new_vect=matmul(matx,vect)
    end function 
    function ftest(xin) result(y)
        implicit none
        real*8, intent(in)::  xin(:,:)
        real*8 ,allocatable :: y(:)
        allocate(y(size(xin,1)))
        y=xin(:,1)**2+xin(:,2)**2
    end function ftest
    subroutine kron_d_cols(A,B,ma,na,nb,AB)

        ! use omp_lib
        implicit none
        
        
        INTEGER, INTENT(IN) :: na,nb,ma
        real*8, DIMENSION(ma,na), INTENT(IN) :: A
        real*8, DIMENSION(ma,nb), INTENT(IN) :: B
        ! real*8, DIMENSION(mb,nb)             :: newB
        real*8, DIMENSION(ma,na*nb), INTENT(out) :: AB
        INTEGER,save :: jj,kk
        ! integer :: numtr
        !  numtr=omp_get_max_threads();
        !  open(10,file='maxthreads.dat')
        !  write(10,*) numtr
        !  close(10)
        AB=0.0d0;
        ! !$OMP do
        ! ! do kk=1,ma*mb
        ! ! AB(kk,kk)=1.0
        ! ! enddo
        ! !$OMP end do
        
        !$OMP parallel do shared(A,B,AB,na,ma)
         do kk=1,na
        
             do jj=1,nb
        
                !call mulscalmx(A(jj,kk),B,mb,nb,newB)
                AB(:,jj+nb*(kk-1))=A(:,kk)*B(:,jj);!newB;       !(A(jj,kk)*B)
                end do
        
        end do
        !$OMP END parallel DO
        
        return
        end subroutine
        SUBROUTINE m_polynomial(xin,yout)!( ngrid_+1, ord, D_st,B1)
            ! CONSTRUCT MONOMIAL BASES
            implicit none 
            real*8, intent(in) :: xin(0:ngrid_-1)
            real*8, intent(out) :: yout(0:ngrid_-1,0:ord-1)
            integer :: cnt

            yout(:,0)=1.0d0
            do cnt=1,ord-1
                yout(:,cnt)=xin**cnt
            end do 

        END SUBROUTINE
        SUBROUTINE SVD_gl(A,U,S,V,M,Nsz)
            ! "economy size SVD"
            implicit none

            real*8 A(M,Nsz),U(M,M),VT(Nsz,Nsz),S(Nsz),V(Nsz,Nsz)
      !
      ! Program computes the matrix singular value decomposition. 
      ! Using Lapack library.
      !
      ! Programmed by sukhbinder Singh
      ! 14th January 2011
      !      
      
      
            real*8,ALLOCATABLE :: WORK(:)
            INTEGER LDA,M,Nsz,LWORK,LDVT,INFO,j,I,LDU
            CHARACTER  JOBU, JOBVT
            real*8 :: iwork(10*Nsz)
            JOBU='S'
            JOBVT='A'
            LDA=M
            LDU=M
            LDVT=Nsz
            
            lwork=-1
            allocate(work(0:0))
            call dgesdd(JOBU,  M, Nsz, A, LDA, S, U, LDU, VT, LDVT,WORK, LWORK,iwork, INFO )
            
            ! LWORK=MAX(1,3*MIN(M,Nsz)+MAX(M,Nsz),5*MIN(M,Nsz))
            lwork=int(work(0))

            deallocate(work)

            ALLOCATE(work(0:lwork-1))
            call dgesdd(JOBU,  M, Nsz, A, LDA, S, U, LDU, VT, LDVT,WORK, LWORK,iwork, INFO )
            ! (jobu, msz, nnsz, VV, msz, Ssvd, Usvd, ldu, Vsvdt,ldvt, work, lwork, iwork, info)
            ! CALL DGESVD(JOBU, JOBVT, M, Nsz, A, LDA, S, U, LDU, VT, LDVT,WORK, LWORK, INFO )
            deallocate(work)        
            V=transpose(VT)
            END SUBROUTINE
end module globals_4test_mod 