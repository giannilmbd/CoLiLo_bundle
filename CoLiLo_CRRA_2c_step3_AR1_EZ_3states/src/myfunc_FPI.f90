subroutine myfunc_FPI(Numinp, x0, fval, parin0)

    use my_toolbox
    use globals

    implicit none
    integer, intent(in):: Numinp
    real(kind=wp), intent(inout) :: x0(Numinp)
    real(kind=wp), intent(out) :: fval(Numinp)
    real(kind=wp), intent(in) :: parin0(7)
    real(kind=wp) :: parin01(4), exit_crit, resid_change
    real(kind=wp) :: X1plus, X2plus
    real(kind=wp) ::  v1plus, v2plus, v1plusm1, v2plusm1,s1plus,A1plus1,A2plus1
    real(kind=wp) ::  xin_V1, xin_V2,Xin_S, bis_a, bis_b, resid_a, resid_b!,xin_S
    real(kind=wp) :: X_n1, X_n2, kappa, C1, C2, C1plus, mu, YWorld, pf, pf_new, L1, L2, ph, Q, Y1, res, pred
    integer :: is_p1, is_p2, ik1_n,ik2_n,ik3_n, iA1_n, iA2_n, is1_n, is2_n, ik1, iter_2, cnt
    real(kind=wp) :: temp_arr1(NS), temp_arr2(NS), possible_pf(nkappas + 1)
    real(kind=wp) :: spline_array1(NS, NS), spline_array2(NS, NS), v1_arr(NS, NS), v2_arr(NS, NS)


    ! control variables (policies)

    ! xin_S=x_input(3)
    ! relative MUC
    ! State variables

    is1_n = INT(parin0(3)) ! current shock index
    is2_n = INT(parin0(4)) ! current shock index
    ik2_n = INT(parin0(5))
    ik3_n = INT(parin0(6))
    ik1_n = INT(parin0(7))

    kappa = kappaS(ik1_n)



   fval(:) = 1.0e6_wp
    iter_2 = 0
    res = 10e6_wp
    ! ITERATE TILL CONVERGENCE
    pf = 1.0d0
    do while (((fval(1) > sig) .or. (fval(2) > sig) .or. (fval(3) > sig)) .and. (iter_2 < itermax))
        iter_2 = iter_2 + 1
        xin_V1 = (x0(1))
        xin_V2 = (x0(2))
        xin_S =  (x0(3))
       !!! BISECTION
        possible_pf = 1.0e0_wp
        ! note that pf must be pf>(1-nu)**(1/(trade_elast-1)) (if trade_elast>1, else it's not so important if nu<<1): see ph expression
        possible_pf(1)=price_down;
        possible_pf(2)=price_up
        bis_a = possible_pf(1)
        bis_b = possible_pf(2)
        cnt = 1
        exit_crit = 100.0e0_wp
        resid_change = exit_crit
A1plus1=A1(ik2_n)**rho_A1*exp(eta1(is1_n))
A2plus1=A2(ik3_n)**rho_A2*exp(eta2(is2_n))

        if (trade_elast < 100.0e0_wp) then
          bisection: do while ((cnt .le. nkappas) .and. (abs(exit_crit) > sig_pf) .and. (abs(resid_change)/abs(exit_crit) > sig_pf))
                if (cnt <= 2) then
                    pf = possible_pf(cnt)
                else
                    pf = (bis_a + bis_b)/2.0e0_wp
                    possible_pf(cnt) = pf
                end if
                ! write (*, '(A55,2x,f10.5)') '[32m pf Value Inputed: [0m ', pf
                call core_model(pf=pf, kappa=kappa, C=C1, Cs=C2, L=L1, Ls=L2, &
                                ph=ph, Q=Q, Y1=Y1, res=res, pred=pred, D=A1plus1, Ds=A2plus1)
                resid_change = exit_crit - res
                exit_crit = res

                if ((abs(res) < sig_pf) .or. (cnt == nkappas)) then
                    exit

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
            end do bisection
        else
            pf = 1.0e0_wp
            call core_model(pf=pf, kappa=kappa,  C=C1, Cs=C2, L=L1, Ls=L2, &
                            ph=ph, Q=Q, Y1=Y1, res=res, pred=pred, D=A1plus1, Ds=A2plus1)
        end if

        !  print*,'stop'
! read(*,*)cnt
        if (cnt >= nkappas) then
            print *, '\033[91m PRICE NOT SOLVED \033[0m'
            print *, 'residual', res
        end if
        !########################################################################
        !  print*,'pf',pf,'redidual',res

        ! call core_model(pf=pf, kappa=kappa, C=C1, Cs=C2, L=L1, Ls=L2,&
        ! D=parin01(2),Ds=parin01(3))
        ! future expectations depend on the update of the policy function Xin_V1 and xin_V2 via OmegaSplus
        Ev1(ik1_n,ik2_n,ik3_n, is1_n, is2_n) = 0e0_wp
        Ev2(ik1_n,ik2_n,ik3_n, is1_n, is2_n) = 0e0_wp
        Es1(ik1_n,ik2_n,ik3_n, is1_n, is2_n) = 0e0_wp
        !!!$omp parallel shared(Ev1,Ev2) private(v1plus,v2plus) default(shared)
   !! !$omp do reduction(+:Ev1,Ev2)
        do is_p1 = 1, NS
            do is_p2 = 1, NS
                v1plus = spline_eval([min(kappa,kappa_u),A1(ik2_n),A2(ik3_n)], &
                                     coeff_v1(:,:,:, is_p1, is_p2), LB_all, UB_all)
                v2plus = spline_eval([min(kappa,kappa_u),A1(ik2_n),A2(ik3_n)], &
                                     coeff_v2(:, :,:, is_p1, is_p2), LB_all, UB_all)
                s1plus= spline_eval([min(kappa,kappa_u),A1(ik2_n),A2(ik3_n)], &
                                    coeff_s1(:,:,:, is_p1, is_p2), LB_all, UB_all)                                     

                ! extrapolation if at the bound(s)

                ! extrapolation if at the bound(s)
                ! ! print*,"line698 test.f90 must allow for undervalue and for over-under-value shocks"
                ! if ((kappa > kappa_u)) then
                !     ! print*,'here3'
                !     v1plus = v1plus + (v1(NK1, is_p1, is_p2) - v1(NK1 - 1,  is_p1, is_p2)) &
                !              /(kappaS(NK1) - kappaS(NK1 - 1))*(kappa - kappa_u)
                !     v2plus = v2plus + (v2(NK1, is_p1, is_p2) - v2(NK1 - 1,  is_p1, is_p2)) &
                !     /(kappaS(NK1) - kappaS(NK1 - 1))*(kappa - kappa_u)
                !     s1plus=s1plus+(s1(NK1, is_p1, is_p2) - s1(NK1 - 1,  is_p1, is_p2)) &
                !     /(kappaS(NK1) - kappaS(NK1 - 1))*(kappa - kappa_u)

                ! end if
                Ev1(ik1_n,ik2_n,ik3_n,  is1_n, is2_n) = Ev1(ik1_n,ik2_n,ik3_n,  is1_n, is2_n) + pi1(is1_n, is_p1)*pi2(is2_n, is_p2)*v1plus ! See Rudebush and Swansson on changing sign with 1/gamma>1
                Ev2(ik1_n,ik2_n,ik3_n,  is1_n, is2_n) = Ev2(ik1_n,ik2_n,ik3_n,  is1_n, is2_n)+ pi1(is1_n, is_p1)*pi2(is2_n, is_p2)*v2plus
                Es1(ik1_n,ik2_n,ik3_n,  is1_n, is2_n) = Es1(ik1_n,ik2_n,ik3_n,  is1_n, is2_n)+ pi1(is1_n, is_p1)*pi2(is2_n, is_p2)*s1plus

            end do
        end do
   !! !$omp END DO
    !!!$omp END PARALLEL
        ! get first order condition
        x0(1) = ((1e0_wp - beta)*(C1**(1e0_wp - 1e0_wp/gamma)/(1e0_wp - 1e0_wp/gamma)-chhi*L1**(1.0e0_wp+varphi)/(1.0e0_wp+varphi))&
         + beta*Ev1(ik1_n,ik2_n,ik3_n, is1_n, is2_n)) ! check sidn Rud&Swans
        x0(2) = ((1e0_wp - beta)*(C2**(1e0_wp - 1e0_wp/gamma)/(1e0_wp - 1e0_wp/gamma)-chhi*L2**(1.0e0_wp+varphi)/(1.0e0_wp+varphi))&
         + beta*Ev2(ik1_n,ik2_n,ik3_n, is1_n, is2_n))
        X0(3) = C1**(-1.0_wp/gamma)*(C1-ph*Y1)+beta*Es1(ik1_n,ik2_n,ik3_n,  is1_n, is2_n)
! write(*,'(a)') ANSI_GREEN//'******************************************************************************************'//ANSI_RESET
!  print*,'CONSUMPTION',C1
!  print*,'LABOR',L1
!  write(*,'(a, 3e20.10)'),'RHS',X0
! !  read(*,*) cnt
!  write(*,'(a)') ANSI_GREEN//'******************************************************************************************'//ANSI_RESET

        fval(1) = abs((x0(1) - xin_V1)/max(xin_V1,1.0e-10_wp))
        fval(2) = abs((x0(2) - xin_V2)/max(xin_V2,1.0e-10_wp))
        fval(3) = abs((X0(3) - Xin_S))

        if ((fval(1) .le. sig) .and. (fval(2) .le. sig).and. (fval(3) .le. sig)) then
            C1v(ik1_n,ik2_n,ik3_n,  is1_n, is2_n) = C1
            L1v(ik1_n,ik2_n,ik3_n,  is1_n, is2_n) = L1
            Y1v(ik1_n,ik2_n,ik3_n,  is1_n, is2_n) = Y1
            pf_sols(ik1_n,ik2_n,ik3_n, is1_n, is2_n) = pred
            phv(ik1_n,ik2_n,ik3_n, is1_n, is2_n) = ph
            RETURN
        else
            x0(1) = smooth*x0(1) + (1 - smooth)*xin_V1
            x0(2) = smooth*x0(2) + (1 - smooth)*xin_V2
            x0(3) = smooth*x0(3) + (1 - smooth)*Xin_S
            ! pf=smooth_pf*pred+(1-smooth_pf)*pf
            ! print*,'vs',x0,'pf',pf
        end if
    end do
! WE ARE GOING TO USE THE FOLLOWING IN ARROW SECURITIES
    C1v(ik1_n,ik2_n,ik3_n,  is1_n, is2_n) = C1
    L1v(ik1_n,ik2_n,ik3_n,  is1_n, is2_n) = L1
    Y1v(ik1_n,ik2_n,ik3_n, is1_n, is2_n) = Y1
    pf_sols(ik1_n,ik2_n,ik3_n, is1_n, is2_n) = pf

    phv(ik1_n,ik2_n,ik3_n,  is1_n, is2_n) = ph
   

end subroutine
