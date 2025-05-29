program test_module

use writeformat, ONLY : wrtfmt

implicit none
!   external  wrtfmt
character*30 :: st,digs,decs
integer :: m,dig,dec
double precision , dimension(10,10) :: A
A=1;
dig=10
dec=10;
m=10;

write(unit=digs,fmt='(i5)') dig
write(unit=decs,fmt='(i5)') dec
! write(*,*) '(f',trim(adjustl(digs)),'.',trim(adjustl(decs)),')'
st=wrtfmt(m,10,5);
! write(*,fmt='(27A)') '(f',trim(adjustl(digs)),'.',trim(adjustl(decs)),')'
write(*,*) st
write(*,trim(st)) A
end program