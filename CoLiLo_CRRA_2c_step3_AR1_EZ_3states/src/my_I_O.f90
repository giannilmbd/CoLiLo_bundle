module my_I_O

contains
    pure function StrNumStr(dir, str1, cnt, str2)

        implicit none
        character(len=*), intent(in) :: str1, str2
        character(len=*), intent(in), optional :: dir
        character(len=100) ::  StrNumStr
        character(len=100) :: namefile1, dir2
        integer, intent(in) :: cnt
        integer :: pos

        write (namefile1, '(i5)') cnt
        namefile1 = adjustr(namefile1)
        pos = index(namefile1, ' ', BACK=.true.)
        if (present(dir) .neqv. .true.) then
            StrNumStr = trim(str1)//trim(namefile1(pos + 1:len(namefile1)))//trim(str2)
        else
            StrNumStr = trim(dir2)//'//'//trim(str1)//trim(namefile1(pos + 1:len(namefile1)))//trim(str2)
        end if

! StrNumStr=trim(str1)//trim(namefile1)//trim(str2)

    end function

    function StrStrStr(dir, str1, cnt, str2) ! just to generalize same formats, but pretty redundant
        ! USE IFPORT
        implicit none
        character(len=*), intent(in) :: str1, str2, cnt
        character(len=*), intent(in), optional :: dir
        character(len=100) ::  StrStrStr
        character(len=100) :: dir2

! integer :: pos

! write(namefile1,'(i5)')cnt
! namefile1=adjustr(namefile1)
! pos=index(namefile1,' ',BACK=.true.)
        if (present(dir) .neqv. .true.) then
            StrStrStr = trim(str1)//trim(cnt)//trim(str2)
        else
                
            dir2=trim(dir)

            StrStrStr = trim(dir2)//'//'//trim(str1)//trim(cnt)//trim(str2)
        end if

! StrNumStr=trim(str1)//trim(namefile1)//trim(str2)

        end function

    end module
