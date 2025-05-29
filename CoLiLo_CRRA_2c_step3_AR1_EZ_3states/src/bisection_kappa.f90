subroutine bisection_kappa(kappa,res)
use my_toolbox, only: spline_eval
use my_kinds_mod, only: wp
use globals
implicit none 

integer:: cnt
real(wp),intent(out) :: kappa,res
real(wp) :: exit_crit,resid_change,possible_k(nkappas)

real(wp) :: bis_a,bis_b,resid_a,resid_b,arrow

bis_a=kappaS(0)
bis_b=kappaS(NK1)
cnt=1
exit_crit=1.0e10
resid_change=1.0e10
      bisection: do while ((cnt .le. nkappas) .and. (abs(exit_crit) > sig_pf) .and. (abs(resid_change)/abs(exit_crit) > sig_pf))
         if (cnt <= 2) then
            kappa = possible_k(cnt)
         else
            kappa = (bis_a + bis_b)/2.0e0_wp
            possible_k(cnt) = kappa
         end if

         ! Compute residual at the current midpoint
arrow = spline_eval([kappa,A1(pos_zeroA1),A2(pos_zeroA2),OmegaS(pos_zeroOm)],coeff_s1(:,:,:,:,pos_zero1,pos_zero2))         
         res=arrow
         resid_change = exit_crit - res
         exit_crit = res

         if ((abs(res) < sig_pf) .or. (cnt == nkappas)) then
            exit
         end if

         ! Update bounds based on the sign of the residual
         if (cnt == 1) then
            resid_a = res
         elseif (cnt == 2) then
            resid_b = res
            if (verbose.and.(sign(1.0e0_wp, resid_a) .eq. sign(1.0e0_wp, resid_b))) then
               print *, "WARNING: Residuals at bounds still have the same sign. Expanding bounds."
               bis_b = bis_b + step_kappa
            end if
         else
            if (sign(1.0e0_wp, resid_a) .eq. sign(1.0e0_wp, res)) then
               resid_a = res
               bis_a = kappa
            elseif (sign(1.0e0_wp, resid_b) .eq. sign(1.0e0_wp, res)) then
               resid_b = res
               bis_b = kappa
            end if
         end if
         cnt = cnt + 1
      end do bisection
   end subroutine bisection_kappa