program run_this
    ! this program runs the cmake command
    implicit none

    call execute_command_line('sh run_all.sh 2&> session.log & ')

end program run_this

! compile this with simply gfortran run_this.f90 -o run_this 