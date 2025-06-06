MODULE SF_VERSION
  implicit none

  !GIT VERSION
  character(len=41),parameter,public :: scifor_version_sha1 = "28b6ed8da5b6117adac1d5c1c6e0b32a468f115b"

contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : print actual version of the software (if any)
  !+-------------------------------------------------------------------+
  subroutine scifor_version()
    integer(4),dimension(8)                  :: dummy
    integer(4)                               :: year
    integer(4)                               :: mese
    integer(4)                               :: day
    integer(4)                               :: h
    integer(4)                               :: m
    integer(4)                               :: s
    integer(4)                               :: ms
    character(len=9),parameter,dimension(12) :: month = (/ &
         'January  ', 'February ', 'March    ', 'April    ', &
         'May      ', 'June     ', 'July     ', 'August   ', &
         'September', 'October  ', 'November ', 'December ' /)
    write(*,"(A)")("SCIFOR VERSION (GIT): "//trim(adjustl(trim(scifor_version_sha1))))
    call date_and_time(values=dummy)
    year = dummy(1)
    mese = dummy(2)
    day  = dummy(3)
    h    = dummy(5)
    m    = dummy(6)
    s    = dummy(7)
    ms   = dummy(8)
    write(*,"(A,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3)")&
         "Timestamp: +",day,trim(month(mese)),year, h,':',m,':',s,'.',ms
    write(*,*)""
    open(10,file="scifor_version.inc")
    write(10,"(A)")"SCIFOR VERSION (GIT): "//trim(adjustl(trim(scifor_version_sha1)))
    write(10,"(A,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3)")&
         "Timestamp: +",day,trim(month(mese)),year, h,':',m,':',s,'.',ms
    write(10,*)""
    close(10)
  end subroutine scifor_version


  subroutine code_version(revision,name)
    character(len=*)                         :: revision
    character(len=*),optional                :: name
    character(len=100)                       :: name_
    integer(4),dimension(8)                  :: dummy
    integer(4)                               :: year
    integer(4)                               :: mese
    integer(4)                               :: day
    integer(4)                               :: h
    integer(4)                               :: m
    integer(4)                               :: s
    integer(4)                               :: ms
    character(len=9),parameter,dimension(12) :: month = (/ &
         'January  ', 'February ', 'March    ', 'April    ', &
         'May      ', 'June     ', 'July     ', 'August   ', &
         'September', 'October  ', 'November ', 'December ' /)
    name_="code";if(present(name))name_=(trim(adjustl(trim(name))))
    write(*,"(A)")(trim(name_)//" VERSION: "//trim(adjustl(trim(revision))))
    call date_and_time(values=dummy)
    year = dummy(1)
    mese = dummy(2)
    day  = dummy(3)
    h    = dummy(5)
    m    = dummy(6)
    s    = dummy(7)
    ms   = dummy(8)
    write(*,"(A,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3)")&
         "Timestamp: +",day,trim(month(mese)),year, h,':',m,':',s,'.',ms
    write(*,*)""
    open(10,file=trim(name_)//"_version.inc")
    write(10,"(A)")trim(name_)//" VERSION: "//trim(adjustl(trim(revision)))
    write(10,"(A,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3)")&
         "Timestamp: +",day,trim(month(mese)),year, h,':',m,':',s,'.',ms
    write(10,*)""
    close(10)
  end subroutine code_version

END MODULE SF_VERSION
