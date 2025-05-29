module find_kappa_mod2
    use my_kinds_mod, only: wp
    use globals,only: pos_zeroA1,pos_zeroA2,pos_zeroOm,pos_zero1,pos_zero2,ANSI_RESET,ANSI_RED,ANSI_GREEN,ANSI_YELLOW,sig_kappa,kappa_min,kappa_max,s1
    implicit none
    contains
    subroutine find_zero_arrow_bisection(xkappa, nres)

        implicit none

        real(kind=wp), intent(out) :: xkappa, nres
        real(kind=wp) :: epsilon,min_val,max_val
        integer :: max_iter, iter
        real(kind=wp) :: a, b, fa, fb, xm, fm
        real(kind=wp) :: step_kappa
    
        ! Set convergence criteria
        epsilon = sig_kappa
        max_iter = 1000000
        step_kappa = 0.01_wp  ! Increment for expanding the range
    
        ! Initialize boundaries
        a = kappa_min
        b = kappa_max
        fa = kappafunc_arrow(a)
        fb = kappafunc_arrow(b)
    
        ! Check if the initial interval contains a zero
        if (sign(1.0e0_wp, fa) .eq. sign(1.0e0_wp, fb)) then
            write(*,*) ANSI_YELLOW//'WARNING: Initial interval does not contain a zero. Proceeding with iterations.'//ANSI_RESET
           
        end if
    
        ! Bisection loop
        do iter = 1, max_iter
            ! Check if bounds still have the same sign
            if (sign(1.0e0_wp, fa) .eq. sign(1.0e0_wp, fb)) then
                write(*,*) ANSI_YELLOW//'WARNING: Residuals at bounds still have the same sign. Expanding bounds.'//ANSI_RESET
                if(a<b) then 
                    a = a - step_kappa
                    b = b + step_kappa
                else
                    a = a + step_kappa
                    b = b - step_kappa
                end if
                
                fa = kappafunc_arrow(a)
                fb = kappafunc_arrow(b)
                cycle
            end if
    
            ! Compute midpoint and its residual
            xm = (a + b) / 2.0_wp
            fm = kappafunc_arrow(xm)
    
            ! Check for convergence
            if (abs(fm) < epsilon) then
                xkappa = xm
                nres = fm
                return
            end if
    
            ! Update bounds
            if (sign(1.0e0_wp, fm) .eq. sign(1.0e0_wp, fa)) then
                a = xm
                fa = fm
            else
                b = xm
                fb = fm
            end if
            min_val = min(fa,fb)
            max_val = max(fa,fb)
            write(*,*) ANSI_YELLOW//'Iteration: ', iter, ' Kappa: ', xm, 'Min s1: ', min_val, ' Max s1: ', max_val, ' Width: ', abs(b - a),ANSI_RESET
            if (abs(b - a) < epsilon) then
                write(*,*) ANSI_GREEN//'Bisection method converged within maximum iterations.'//ANSI_RESET
                xkappa = xm
                nres = fm
                return
            end if
        end do
    
        ! Report if max iterations reached without convergence
        write(*,*) ANSI_RED//'Bisection method did not converge within maximum iterations.'//ANSI_RESET
        xkappa = xm
        nres = fm
    
    end subroutine find_zero_arrow_bisection
    

! Function to evaluate the residual
real(kind=wp) function kappafunc_arrow(xkappa)

    real(kind=wp), intent(in) :: xkappa
    real(kind=wp) :: fm
    logical ::  converged
    integer, parameter :: size_x = 1

    call solve_policy_function(xkappa, converged)
    
    if (.not. converged) then
        write(*,*) ANSI_RED//'Policy function did not converge at kappa = ', xkappa,ANSI_RESET
        fm = huge(1.0_wp)
    else
        call solve_arrow_securities()
        fm = s1(pos_zeroA1,pos_zeroA2,pos_zeroOm,pos_zero1,pos_zero2)
    end if
    

    kappafunc_arrow = fm

end function kappafunc_arrow

end module