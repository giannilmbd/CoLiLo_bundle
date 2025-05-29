program main
    use quasi_newton_mod
    implicit none
    integer, parameter :: n = 2, maxit = 100
    real(kind=8), parameter :: tol = 1.0d-8
    real(kind=8), dimension(n) :: x
    integer :: i

    ! Initial guess for the solution
    x = [1.0d0, 1.0d0]

    ! Solve the system of non-linear equations using the Quasi-Newton method
    call quasi_newton(F, x, n, maxit, tol)

    ! Print the solution
    write(*,*) "Solution:"
    do i = 1, n
        write(*,*) "x(",i,") = ", x(i)
    end do

end program main

! ifort -O3 -r8 ../../modules/mybroyden.f90 test_broyden.f90 -o boyd.out