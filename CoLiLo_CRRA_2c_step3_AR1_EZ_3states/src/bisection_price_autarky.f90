subroutine bisection_price_autarky(pf, C1, L1,C2,L2,A1plus,A2plus,ph,Q)
    use globals
    implicit none
    integer :: cnt
    real(kind=wp),intent(in) :: A1plus,A2plus
    real(kind=wp), intent(out) :: ph,pf,C1,L1,C2,L2,Q
    real(kind=wp) :: pred,res,Y1
    real(kind=wp) :: exit_crit, resid_change, bis_a, bis_b, resid_a, resid_b, possible_pf(nkappas + 1)
    ! ################ BISECTION FOR PRICES
    cnt = 1
    possible_pf = 1.0e0_wp
    ! note that pf must be pf>(1-nu)**(1/(trade_elast-1)) (if trade_elast>1, else it's not so important if nu<<1): see ph expression
    possible_pf(1) = price_down
    possible_pf(2) = price_up
    bis_a = possible_pf(1)
    bis_b = possible_pf(2)
    exit_crit = 100.0e0_wp
    resid_change = exit_crit
    if (trade_elast < 100.0e0_wp) then
        fpi_price:do while ((cnt .le. nkappas) .and. (abs(exit_crit) > sig_pf)) ! .and. (abs(resid_change)/abs(exit_crit) > sig_pf)
            if (cnt <= 2) then
                pf = possible_pf(cnt)
            else
                pf = (bis_a + bis_b)/2.0e0_wp
                possible_pf(cnt) = pf
            end if
            ! write (*, '(A55,2x,f10.5)') '[32m pf Value Inputed: [0m ', pf
            call core_model_autarky(pf=pf, C=C1, Cs=C2, L=L1, Ls=L2, &
                                    ph=ph, Q=Q, Y1=Y1, res=res, pred=pred, D=A1plus, Ds=A2plus)
            resid_change = exit_crit - res
            exit_crit = res
            !  write (*, '(A35,2x,f10.5)') '[32m pf iteration percent: [0m', dble(ik1)/dble(nkappas)*100.0e0_wp
            !  write (*, '(A35,2x,f10.5)') '[32m pf value: [0m ', pred
            !  write (*, '(A35,2x,f10.5)') '[36m Current pf Residual: [0m', res
            if ((abs(res) < sig_pf) .or. (cnt == nkappas)) then
                exit fpi_price

            end if

            if (cnt == 1) then
                resid_a = res
            elseif (cnt == 2) then
                resid_b = res
                if (sign(1.0e0_wp, resid_a) .eq. sign(1.0e0_wp, resid_b)) then ! maybe better throw an error. Problem: if process is unsupervised
                    bis_b = bis_b + step_kappa

                    pf = (bis_a + bis_b)/2.0e0_wp
                end if
            else

                if (sign(1.0e0_wp, resid_a) .eq. sign(1.0e0_wp, res)) then
                    resid_a = res
                    bis_a = pf
                elseif (sign(1.0e0_wp, resid_b) .eq. sign(1.0e0_wp, res)) then
                    resid_b = res
                    bis_b = pf
                end if
            end if
            cnt = cnt + 1
        end do fpi_price
    else
        pf = 1.0e0_wp
        call core_model_autarky(pf=pf, C=C1, Cs=C2, L=L1, Ls=L2, &
                                ph=ph, Q=Q, Y1=Y1, res=res, pred=pred, D=A1plus, Ds=A2plus)
    end if
!####################
    if ((verbose).and.(abs(res)>sig_pf))write(*,'(a,e20.10)') ANSI_RED // 'RESIDUAL PRICE AUTARKY ABOVE THRESHOLD' // ANSI_RESET, res
end subroutine bisection_price_autarky
