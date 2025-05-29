real(wp) function F_kappa(x_in) 
use globals
implicit none

real(wp), intent(in) :: x_in

call solve_policy_function()
F_kappa=resid_A0(1)
end function