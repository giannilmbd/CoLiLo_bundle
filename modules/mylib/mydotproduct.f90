subroutine mydotproduct(X,Y,N,Z)

implicit none

double precision ,intent(in):: X(N),Y(N)
double precision ,intent(out):: Z
integer ,intent(in):: N
integer :: i
Z=sum((/(X(i)*Y(i),i=1,N)/))
return
end subroutine
