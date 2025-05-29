program testrnorm

implicit none

double precision :: g
double precision, external :: rnorm

g=rnorm()

write(*,*) g
end program