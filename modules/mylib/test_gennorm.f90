program test_gennorm

implicit none

double precision :: epsilon_
real :: av,sd
real,External :: gennor
integer ::  ISEED1,ISEED2,jj
real :: time0,time1
ISEED1=1967
ISEED2=2007
call SETALL(ISEED1,ISEED2)
av=3.;sd=1.;
open(10,file='data.dat');
write(10,*) 
close(10);
open(10,file='data.dat',position='append');
call cpu_time(time0)
do jj=1,10000
epsilon_=gennor(av,sd)

write(10,*) dble(epsilon_)
enddo
close(10);
call cpu_time(time1)

print *, 'time elapsed = ', time1-time0
! write(*,*) epsilon_
end program
