program test
    use globals_4test_mod 
    use mypolynte_mod
    use ogpf
    use linespaced_mod
    
    
    implicit none
character(len=100) :: form,nchar
    integer :: i,j
    real*8,dimension(0:num_coeffs-1,0:ngrid_-1) :: Oinv
    real*8,dimension(0:ngrid_-1,0:ord-1) :: newvalue1, newvalue2
    real*8,dimension(0:ngrid_-1,0:ord**2-1) :: newvalue
    real*8,dimension(1:ngrid_,2):: X4f
    type(params) :: param
    TYPE(gpf):: gp
 ! real*8, allocatable ::work(:)
 ! real*8 :: work0(1),Dbi(0:num_coeffs),Ebi(0:num_coeffs-1)
 ! real*8 :: tauq(0:num_coeffs), taup(0:num_coeffs)

! points on the grid
    D_st(0:ngrid_-1)=linespaced_values(ngrid_,-0.5d0,0.5d0)
    D_stst(0:ngrid_-1)=linespaced_values(ngrid_,-0.5d0,0.5d0)
    Ds_st(0:ngrid_-1)=linespaced_values(ngrid_,-0.5d0,0.5d0)
    Ds_stst(0:ngrid_-1)=linespaced_values(ngrid_,-0.5d0,0.5d0)
    param%D=linespaced_values(ngrid_,-0.5d0,0.5d0)
    param%Ds=linespaced_values(ngrid_,-0.5d0,0.5d0)
    ! Ones(0:2,0:2)=reshape(source=(/01.0d0,0.0d0,0.0d0,1.0d0/),shape=shape(Ones))
    ! B1(0:2,0:2)=reshape(source=(/1.0d0,0.0d0,0.0d0,1.0d0,0.5d0,0.0d0,1.0d0,0.5d0,0.25d0/),shape=shape(B1))
    ! B2(0:2,0:2)=reshape(source=(/1.0d0,0.0d0,0.0d0,1.0d0,0.5d0,0.0d0,1.0d0,0.5d0,0.25d0/),shape=shape(B2))
    
! print function values on the grid
    X4f(:,1)=param%D
    X4f(:,2)=param%Ds
    write(*,*) '\033[1;31m X4f \033[0m \n',X4f (:,1)
    write(*,*) '\033[1;31m Function \033[0m \n',ftest(X4f)
    CALL gp%title('pf \\& XX Residual')
CALL gp%xlabel('points')
CALL gp%ylabel('ftest')
Call gp%options('set style data linespoints')
!Call Plot to draw a vector against a vector of data
CALL gp%plot(X4f(:,1), ftest(x4f))

! obtain the (X*X')^-1*X' matrix
    call myinterpol(Oinv)
! evaluate the interpolatioin
    do i=1,ord

    newvalue1(:,i)=X4f(:,1)**i
    newvalue2(:,i)=X4f(:,2)**i
    enddo 

    !###################### MUST SIMPLY DO COEFFICIENTS * VALUES TO GET A VECTOR
    call kron_d_cols(newvalue1,newvalue2,ngrid_+1,ord,ord,newvalue);
    newfunc=matmul(param%OnGrid,newvalue)

! Convert integer ngrid_ to a character string
    write(nchar, '(I0)') size(Oinv,2)
    nchar = adjustl(trim(nchar))
form='(1x,'//trim(nchar)//'e12.6)'   
! print*,form  
    write(*,trim(form)) Oinv
   

end program test



!ifort -o test ../chebyshev/chebyshev_polynomial.f90 mylib/solve_lapack.f90 test.f90 mypolynte_mod.f90 globals_4test_mod.f90 -L/usr/lib64/ -llapack -lblas