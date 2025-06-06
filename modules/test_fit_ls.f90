program test_fit_ls
use my_kinds_mod, only:wp
implicit none

integer, parameter:: n=5,m=2
real(wp) :: x(n,m),y(n), std_err,beta(m+1)

call random_number(x)
call random_number(y)

call fit_ls(x,y,beta,std_err)

print*, 'beta',beta
print*, 'std_err',std_err

end program 

! ifort ../my_kinds_mod.f90 fit_ls.f90 test_fit_ls.f90 -llapack -o test_fit_ls.out

