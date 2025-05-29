module find_zero_bisection_mod
    use toolbox, only: spline_eval
    use my_kinds_mod, only: wp
    use globals,only: coeff_s1,pos_zero1,pos_zero2,LB_all,UB_all,ANSI_RESET,ANSI_RED,sig_kappa
    implicit none
    contains
subroutine find_zero_bisection(xkappa_min, xkappa_max, xkappa, nres)

    implicit none
    real(kind=wp), intent(in) :: xkappa_min, xkappa_max
    real(kind=wp), intent(out) :: xkappa, nres
    real(kind=wp) :: epsilon
    integer :: max_iter, iter
    real(kind=wp) :: a, b, fa, fb, xm, fm

    ! Set convergence criteria
    epsilon = sig_kappa
    max_iter = 10000

    ! Initialize boundaries
    a = xkappa_min
    b = xkappa_max
    fa = kappafunc(a)
    fb = kappafunc(b)

    ! Check if the initial interval contains a zero
    if (sign(1.0e0_wp, fa) .eq. sign(1.0e0_wp, fb)) then
        write(*,*) ANSI_RED//'Initial interval does not contain a zero.'//ANSI_RESET
        return
    end if

    ! Bisection loop
    do iter = 1, max_iter
        xm = (a + b) / 2.0_wp
        fm = kappafunc(xm)

        if (abs(fm) < epsilon) then
            xkappa = xm
            nres = fm
            return
        end if

        if (sign(1.0e0_wp, fm) .eq. sign(1.0e0_wp, fa)) then
            a = xm
            fa = fm
        else
            b = xm
            fb = fm
        end if
    end do

    ! Report if max iterations reached without convergence
    write(*,*) ANSI_RED//'Bisection method did not converge within maximum iterations.'//ANSI_RESET
    xkappa = xm
    nres = fm

end subroutine find_zero_bisection

! Function to evaluate the residual
real(kind=wp) function kappafunc(xkappa)

    real(kind=wp), intent(in) :: xkappa
    real(kind=wp) :: results(1), nres
    integer, parameter :: size_x = 1

    
    results = spline_eval(xkappa, coeff_s1(:, pos_zero1, pos_zero2), LB_all, UB_all)

    kappafunc = results(1)
end function kappafunc
end module