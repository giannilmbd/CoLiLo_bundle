module writeformat
 implicit none
save
! interface
contains
 function wrtfmt(pm,dig,dec)
character*27:: wrtfmt
  character*5 :: sz
  character*5 :: digs,decs
integer,intent(in) :: pm,dig,dec

write(unit=digs,fmt='(i5)') dig
write(unit=decs,fmt='(i5)') dec
write(Unit=sz,fmt='(i5.0)') pm

! write(*,*) len(sz),sz
write(unit=wrtfmt,fmt='(27A)') '(',sz,'f',trim(adjustl(digs)),'.',trim(adjustl(decs)),')'
wrtfmt=trim(wrtfmt)
! write(*,*) trim(szl),len(szl)
!  write(*,trim(szl)) A
end function wrtfmt

! end interface
end module writeformat

!  ifort -c mexf90.f90