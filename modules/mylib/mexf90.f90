module mexf90
save
interface

  function mxGetPr(pm)
  integer,pointer :: mxGetPr
integer :: pm
end function mxGetPr
function mxGetM(pm)
integer :: mxGetM,pm
end function mxGetM
function mxGetN(pm)
integer :: mxGetN,pm
end function mxGetN
function mxCreateDoubleMatrix(m,n,type)
integer :: m,n,type,mxCreateDoubleMatrix
end function mxCreateDoubleMatrix
function mxGetScalar(pm)
integer :: pm
double precision :: mxGetScalar
end function mxGetScalar
end interface
end module mexf90

!  ifort -c mexf90.f90