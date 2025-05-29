module find_kappa_mod
    use my_kinds_mod, only: wp
    use globals, only: sig_kappa, ANSI_RESET, ANSI_RED, ANSI_YELLOW, ANSI_GREEN,kappa_min,kappa_max, s1, pos_zeroA1, pos_zeroA2, pos_zeroOm, pos_zero1, pos_zero2
    implicit none

    private
    public :: find_optimal_kappa

contains
    subroutine find_optimal_kappa(kappa_opt, final_s1)
        implicit none
        real(kind=wp), intent(out) :: kappa_opt, final_s1
        real(kind=wp) :: epsilon
        integer :: max_iter, iter, expand_count
        real(kind=wp) :: a, b, fa, fb, xm, fm
        real(kind=wp) :: step_kappa
        logical :: converged
        integer, parameter :: max_expansions = 50
    
        ! Set convergence criteria
        epsilon = sig_kappa
        max_iter = 1000000
        step_kappa = 0.1_wp  ! Initial step size for expanding bounds
        expand_count = 0
    
        ! Initialize boundaries
        a = kappa_min
        b = kappa_max
        
        ! First solve the model and get s1 values at boundaries
        call solve_policy_function(a, converged)
        if (.not. converged) then
            write(*,*) ANSI_RED//'Policy function did not converge at kappa = ', a,ANSI_RESET
            fa = huge(1.0_wp)
        else
            call solve_arrow_securities(a)
            fa = s1(pos_zeroA1,pos_zeroA2,pos_zeroOm,pos_zero1,pos_zero2)
        end if

        call solve_policy_function(b, converged)
        if (.not. converged) then
            write(*,*) ANSI_RED//'Policy function did not converge at kappa = ', b,ANSI_RESET
            fb = huge(1.0_wp)
        else
            call solve_arrow_securities(b)
            fb = s1(pos_zeroA1,pos_zeroA2,pos_zeroOm,pos_zero1,pos_zero2)
        end if
    
        ! Bisection loop
        iteration: do iter = 1, max_iter
            ! Check if bounds still have the same sign
            if (sign(1.0_wp, fa) .eq. sign(1.0_wp, fb)) then
                expand_count = expand_count + 1
                if (expand_count > max_expansions) then
                    write(*,*) ANSI_RED//'Failed to bracket root after ', max_expansions, ' expansions.'//ANSI_RESET
                    write(*,*) ANSI_RED//'Current bounds: [', a, ', ', b, '] with values [', fa, ', ', fb, ']'//ANSI_RESET
                    exit
                end if

                write(*,*) ANSI_YELLOW//'WARNING: s1 values at bounds have the same sign. Expanding bounds.'//ANSI_RESET
                ! Expand bounds outward from their current positions
                a = a - sign(1.0_wp, fa) * step_kappa  ! Move a opposite to fa's sign
                b = b + sign(1.0_wp, fb) * step_kappa  ! Move b in fb's sign direction
                step_kappa = step_kappa * 2.0_wp  ! Double the step size for next expansion
                
                call solve_policy_function(a, converged)
                if (.not. converged) then
                    write(*,*) ANSI_RED//'Policy function did not converge at kappa = ', a,ANSI_RESET
                    fa = huge(1.0_wp)
                else
                    call solve_arrow_securities(a)
                    fa = s1(pos_zeroA1,pos_zeroA2,pos_zeroOm,pos_zero1,pos_zero2)
                end if

                call solve_policy_function(b, converged)
                if (.not. converged) then
                    write(*,*) ANSI_RED//'Policy function did not converge at kappa = ', b,ANSI_RESET
                    fb = huge(1.0_wp)
                else
                    call solve_arrow_securities(b)
                    fb = s1(pos_zeroA1,pos_zeroA2,pos_zeroOm,pos_zero1,pos_zero2)
                end if
                
                write(*,'(A,2E15.6)') 'New bounds: ', a, b
                write(*,'(A,2E15.6)') 'New values: ', fa, fb
                cycle
            end if
    
            ! Compute midpoint
            xm = (a + b) / 2.0_wp
            
            ! Solve model at midpoint
            call solve_policy_function(xm, converged)
            if (.not. converged) then
                write(*,*) ANSI_RED//'Policy function did not converge at kappa = ', xm,ANSI_RESET
                fm = huge(1.0_wp)
            else
                call solve_arrow_securities(xm)
                fm = s1(pos_zeroA1,pos_zeroA2,pos_zeroOm,pos_zero1,pos_zero2)
            end if
    
            ! Check convergence based on both function value and interval width
            if (abs(fm) < epsilon .or. abs(b - a) < epsilon) then
                kappa_opt = xm
                final_s1 = fm
                write(*,*) ANSI_GREEN//'Found optimal kappa = ', kappa_opt, ' with s1 = ', final_s1//ANSI_RESET
                write(*,*) ANSI_GREEN//'Converged after ', iter, ' iterations with interval width = ', abs(b - a)//ANSI_RESET
                return
            end if
    
            ! Update bounds
            if (sign(1.0_wp, fm) .eq. sign(1.0_wp, fa)) then
                a = xm
                fa = fm
            else
                b = xm
                fb = fm
            end if

            ! Print progress every 100 iterations
            if (mod(iter, 100) == 0) then
                write(*,'(A,I6,A,F12.6,A,E12.4,A,E12.4)') 'Iteration: ', iter, &
                    ' Kappa: ', xm, ' s1: ', fm, ' Width: ', abs(b - a)
            end if
        end do iteration
    
        ! Report if max iterations reached without convergence
        write(*,*) ANSI_RED//'Bisection method did not converge within maximum iterations.'//ANSI_RESET
        write(*,*) ANSI_RED//'Final interval width = ', abs(b - a), ' with s1 = ', fm//ANSI_RESET
        kappa_opt = xm
        final_s1 = fm
    
    end subroutine find_optimal_kappa

    
end module find_kappa_mod 