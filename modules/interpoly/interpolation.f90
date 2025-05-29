! this program shows the use of the fgsl library for the polynomial interpolation of a three dimensional funciton
! the program is based on the example of the fgsl library
! the function to be interpolated is the function of the example of the fgsl library
! the function is defined in the subroutine f
program interpolation
use fgsl
implicit none
real(fgsl_double), dimension(3,3) :: x ! the points to be interpolated
real(fgsl_double), dimension(3) :: y ! the values of the function to be interpolated
real(fgsl_double), dimension(3) :: c  ! the coefficients of the polynomial
real(fgsl_double) :: z ! the value of the interpolated function
integer(fgsl_int) :: i, status   ! loop variable
! determine the points to be interpolated
do i=1,3
x(i,1)=i
x(i,2)=i
x(i,3)=i
end do

! determine the values of the function to be interpolated
do i=1,3
call f(x(i,1),x(i,2),x(i,3),y(i))
end do
!interpolation
! determine the coefficients of the polynomial
status= fgsl_poly_dd_int(c,x,y)

! determine the value of the interpolated function
call f(1.5,1.5,1.5,z)    ! the value of the function at the point (1.5,1.5,1.5)
print *,z ! print the value of the function at the point (1.5,1.5,1.5)
! determine the value of the interpolated function
z=fgsl_poly_eval(c,1.5_fgsl_double) ! determine the value of the interpolated function
print *,z ! print the value of the interpolated function
end program interpolation
! the function to be interpolated
subroutine f(x,y,z,out)
    use fgsl
implicit none
real(fgsl_double) :: x,y,z,out
out=x**2+y**2+z**2
end subroutine f

! compile the program with the command
! ifort -o interpolation interpolation.f90 -L/usr/local/lib -I/usr/local/include/fgsl -lfgsl



