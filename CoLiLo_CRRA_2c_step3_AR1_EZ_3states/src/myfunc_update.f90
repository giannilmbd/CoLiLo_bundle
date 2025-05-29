subroutine myfunc_update(Numinp, x0, fval, parin0)
! use bspline_mod
use bspline_oo_module

   use my_toolbox, only:spline_eval
   use globals
   use interp_mod, only: quadlinear_eval,bspline4_eval
   implicit none
   
   integer, intent(in):: Numinp
   real(kind=wp), intent(inout) :: x0(Numinp)
   real(kind=wp), intent(out) :: fval(Numinp)
   real(kind=wp), intent(in) :: parin0(100)
   real(kind=wp) :: parin01(4), exit_crit, resid_change,x_states(3)
   real(kind=wp) :: X1plus, X2plus, OmegaSplus, OmegaSplusplus, Omega1, Omega2, Omega1plus, A1plus1, A2plus1,OmegaS_n
   real(kind=wp) ::  v1plus(NS,NS), v2plus(NS,NS), v1plusm1, v2plusm1,s1plus,A1_n,A2_n
   real(kind=wp) ::  xin_V1, xin_V2, bis_a, bis_b, resid_a, resid_b!
   real(kind=wp) :: X_n1, X_n2, kappa, C1, C2, C1plus, mu, YWorld, pf, pf_new, L1, L2, ph, Q, Y1, res, pred
   integer :: is_p1, is_p2, ik2_n, ik3_n, ik4_n, is1_n, is2_n, iter_2, cnt
   real(kind=wp) :: temp_arr1(NS), temp_arr2(NS), possible_pf(nkappas + 1)
   real(kind=wp) :: spline_array1(NS, NS), spline_array2(NS, NS), v1_arr(NS, NS), v2_arr(NS, NS)



   ! control variables (policies)

   ! xin_S=x_input(3)
   ! relative MUC
   ! State variables



   is1_n = INT(parin0(1)) ! current shock index
   is2_n = INT(parin0(2)) ! current shock index
   kappa = parin0(3)
   ik2_n = INT(parin0(4))
   ik3_n = INT(parin0(5))
   ik4_n = INT(parin0(6))

! parin01(1) = parin0(8)
   A1_n = A1(ik2_n)
   A2_n = A2(ik3_n)
   OmegaS_n = OmegaS(ik4_n)

! compute E_{t-1}V_t, which depends on state t-1 and expected innovations
   Ev1m1(ik2_n, ik3_n,ik4_n, is1_n,is2_n) = 0e0_wp
   Ev2m1(ik2_n, ik3_n,ik4_n, is1_n,is2_n) = 0e0_wp

!SHOULD NOT PARALLELIZE THIS AS THE OUTER LAYER IS PARALLELIZED

   x_states= [A1_n, A2_n,OmegaS_n ]
   do is_p2 = 1, NS
      do is_p1 = 1, NS
      !   v1plusm1=abs(bspline4_eval(x_states, tx_v1(:,is_p1,is_p2), ty_v1(:,is_p1,is_p2),&
      !    tz_v1(:,is_p1,is_p2), tq_v1(:,is_p1,is_p2), bcoef4_v1(:,:,:,:,is_p1,is_p2),  .true.))
      !    v2plusm1=abs(bspline4_eval(x_states, tx_v2(:,is_p1,is_p2), ty_v2(:,is_p1,is_p2),&
      !    tz_v2(:,is_p1,is_p2), tq_v2(:,is_p1,is_p2), bcoef4_v2(:,:,:,:,is_p1,is_p2),  .true.))
! lagged expectations
 
            v1plusm1 = abs(spline_eval(x_states, &
               coeff_v1(:, :, :, is_p1, is_p2), LB_all, UB_all))
            v2plusm1 = abs(spline_eval(x_states, &
               coeff_v2(:, :, :, is_p1, is_p2), LB_all, UB_all))


            


         Ev1m1(ik2_n, ik3_n,ik4_n,is1_n,is2_n) = Ev1m1(ik2_n, ik3_n,ik4_n, is1_n,is2_n) &
            + pi1(is1_n, is_p1)*pi2(is2_n, is_p2)*(v1plusm1)**theta ! See Rudebush and Swansson on changing sign with 1/gamma>1
         Ev2m1(ik2_n, ik3_n,ik4_n, is1_n,is2_n) = Ev2m1(ik2_n, ik3_n,ik4_n, is1_n,is2_n) &
            + pi1(is1_n, is_p1)*pi2(is2_n, is_p2)*(v2plusm1)**theta

      end do
   end do


   fval(:) = 1.0e6_wp
   iter_2 = 0
   res = 10e6_wp
   ! ITERATE TILL CONVERGENCE

   ! do while (((fval(1) > sig) .or. (fval(2) > sig)) .and. (iter_2 < itermax))
   ! iter_2 = iter_2 + 1
   xin_V1 = (x0(1))
   xin_V2 = (x0(2))
   
   !!! BISECTION
  
   ! note that pf must be pf>(1-nu)**(1/(trade_elast-1)) (if trade_elast>1, else it's not so important if nu<<1): see ph expression
   possible_pf(1)=price_down;
   possible_pf(2)=price_up
   bis_a = possible_pf(1)
   bis_b = possible_pf(2)
   cnt = 1
   exit_crit = 100.0e0_wp
   resid_change = exit_crit
   A1plus1=(A1(ik2_n))**rho_A1*exp(eta1(is1_n))
   A2plus1=(A2(ik3_n))**rho_A2*exp(eta2(is2_n))
   Omega1 = beta*(Ev1m1(ik2_n, ik3_n,ik4_n, is1_n,is2_n))**(1.0_wp/theta - 1.0_wp)*(xin_V1)**(theta - 1.0_wp) ! See RUdebusch and Swansson for change of sign

   Omega2 = beta*(Ev2m1(ik2_n, ik3_n,ik4_n, is1_n,is2_n))**(1.0_wp/theta - 1.0_wp)*(xin_V2)**(theta - 1.0_wp)

! update state variables
! OmegaS (ie get OmegaS_t, which enters consumption shares and thus defines allocations)
   OmegaSplus = Omega2/Omega1*OmegaS_n
   if (trade_elast < 100.0e0_wp) then
      call core_model(pf=bis_a, kappa=kappa,OmegaS_in=OmegaSplus, C=C1, Cs=C2, L=L1, Ls=L2, &
         ph=ph, Q=Q, Y1=Y1, res=resid_a, pred=pred, D=A1plus1, Ds=A2plus1)
      call core_model(pf=bis_b, kappa=kappa,OmegaS_in=OmegaSplus, C=C1, Cs=C2, L=L1, Ls=L2, &
         ph=ph, Q=Q, Y1=Y1, res=resid_b, pred=pred, D=A1plus1, Ds=A2plus1)

! Check if initial bounds are valid
      do while(sign(1.0e0_wp, resid_a) .eq. sign(1.0e0_wp, resid_b))
         print *, ANSI_RED//"WARNING: Initial bounds for prices do not bracket the root. Adjusting bounds."//ANSI_RESET
         print*, ANSI_BLUE//"resid_a = ", resid_a, " resid_b = ", resid_b,ANSI_RESET
         ! Optionally, terminate or adjust bounds here.
               if(bis_a<bis_b) then   
               bis_a = max(1.0e-6,(bis_a - step_kappa))
               bis_b = min(20.0_wp,bis_b + step_kappa)
               else 
               bis_a = bis_a +step_kappa
               bis_b = min(20.0_wp,bis_b + step_kappa)
               endif 
                     call core_model(pf=bis_a, kappa=kappa,OmegaS_in=OmegaSplus, C=C1, Cs=C2, L=L1, Ls=L2, &
         ph=ph, Q=Q, Y1=Y1, res=resid_a, pred=pred, D=A1plus1, Ds=A2plus1)
      call core_model(pf=bis_b, kappa=kappa,OmegaS_in=OmegaSplus, C=C1, Cs=C2, L=L1, Ls=L2, &
         ph=ph, Q=Q, Y1=Y1, res=resid_b, pred=pred, D=A1plus1, Ds=A2plus1)
         cycle
            end do

      bisection: do while ((cnt .le. nkappas) .and. (abs(exit_crit) > sig_pf) .and. (abs(resid_change)/abs(exit_crit) > sig_pf))
         if (cnt <= 2) then
            pf = possible_pf(cnt)
         else
            pf = (bis_a + bis_b)/2.0e0_wp
            possible_pf(cnt) = pf
         end if

         ! Compute residual at the current midpoint
         call core_model(pf=pf, kappa=kappa,OmegaS_in=OmegaSplus ,C=C1, Cs=C2, L=L1, Ls=L2, &
            ph=ph, Q=Q, Y1=Y1, res=res, pred=pred, D=A1plus1, Ds=A2plus1)
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
            do while (verbose.and.(sign(1.0e0_wp, resid_a) .eq. sign(1.0e0_wp, resid_b))) 
               print *, "WARNING: Residuals at bounds still have the same sign. Expanding bounds."
            !!!!!!!!   
               if(bis_a<bis_b) then   
               bis_a = max(1.0e-6,(bis_a - step_kappa))
               bis_b = min(20.0_wp,bis_b + step_kappa)
               else 
               bis_a = bis_a +step_kappa
               bis_b = min(20.0_wp,bis_b + step_kappa)
               endif 
                     call core_model(pf=bis_a, kappa=kappa,OmegaS_in=OmegaSplus, C=C1, Cs=C2, L=L1, Ls=L2, &
         ph=ph, Q=Q, Y1=Y1, res=resid_a, pred=pred, D=A1plus1, Ds=A2plus1)
      call core_model(pf=bis_b, kappa=kappa,OmegaS_in=OmegaSplus, C=C1, Cs=C2, L=L1, Ls=L2, &
         ph=ph, Q=Q, Y1=Y1, res=resid_b, pred=pred, D=A1plus1, Ds=A2plus1)
         cycle 
            !!!!!!!!!!!
            
            end do
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
      call core_model(pf=pf, kappa=kappa,OmegaS_in=OmegaSplus,  C=C1, Cs=C2, L=L1, Ls=L2, &
         ph=ph, Q=Q, Y1=Y1, res=res, pred=pred, D=A1plus1, Ds=A2plus1)
   end if

   !  print*,'stop'
! read(*,*)cnt
   if ((verbose).and.(abs(res)>sig_pf)) then
! if (abs(res)>sig_pf) then
      write(*,'(a,e20.10)') ANSI_RED // 'RESIDUAL PRICE CM ABOVE THRESHOLD' // ANSI_RESET, res
   endif

   if (cnt>=nkappas) then
      print*,'\033[91m PRICE NOT SOLVED \033[0m'
      print*,'residual',res
   endif
   !########################################################################


   ! call core_model(pf=pf, kappa=kappa, C=C1, Cs=C2, L=L1, Ls=L2,&
   ! D=parin01(2),Ds=parin01(3))
   ! future expectations depend on the update of the policy function Xin_V1 and xin_V2 via OmegaSplus
   Ev1(ik2_n,ik3_n,ik4_n, is1_n, is2_n) = 0e0_wp
   Ev2(ik2_n,ik3_n,ik4_n, is1_n, is2_n) = 0e0_wp
   !!!$omp parallel shared(Ev1,Ev2) private(v1plus,v2plus) default(shared)
   !! !$omp do reduction(+:Ev1,Ev2)

   x_states=[A1plus1,A2plus1,OmegaSplus]
   do is_p2 = 1, NS
      do is_p1 = 1, NS

            v1plus(is_p1,is_p2) = spline_eval(x_states, &
               coeff_v1(:,:,:, is_p1, is_p2), LB_all, UB_all)
               
            v2plus(is_p1,is_p2) = spline_eval(x_states, &
               coeff_v2(:,:,:,  is_p1, is_p2), LB_all, UB_all)


         Ev1(ik2_n,ik3_n,ik4_n, is1_n, is2_n) = Ev1(ik2_n,ik3_n,ik4_n, is1_n, is2_n)   + pi1(is1_n, is_p1)*pi2(is2_n, is_p2)*v1plus(is_p1,is_p2)**theta ! See Rudebush and Swansson on changing sign with 1/gamma>1
         Ev2(ik2_n,ik3_n,ik4_n, is1_n, is2_n) = Ev2(ik2_n,ik3_n,ik4_n, is1_n, is2_n)   + pi1(is1_n, is_p1)*pi2(is2_n, is_p2)*v2plus(is_p1,is_p2)**theta

      end do
   end do

   !! !$omp END DO
   !!!$omp END PARALLEL
   ! get first order condition
   x0(1) = -(1e0_wp - beta)*(C1**(1e0_wp - 1e0_wp/gamma)/(1e0_wp - 1e0_wp/gamma)-chhi*L1**(1.0e0_wp+varphi)/(1.0e0_wp+varphi)) + beta*(Ev1(ik2_n,ik3_n,ik4_n, is1_n, is2_n))**(1.0e0_wp/theta)
   x0(2) = -(1e0_wp - beta)*(C2**(1e0_wp - 1e0_wp/gamma)/(1e0_wp - 1e0_wp/gamma)-chhi*L2**(1.0e0_wp+varphi)/(1.0e0_wp+varphi)) + beta*(Ev2(ik2_n,ik3_n,ik4_n, is1_n, is2_n))**(1.0e0_wp/theta)


   fval(1) = 2.0_wp*abs((x0(1) - xin_V1)/(xin_V1+x0(1)))
   fval(2) = 2.0_wp*abs((x0(2) - xin_V2)/(xin_V2+x0(2)))

! WE ARE GOING TO USE THE FOLLOWING IN ARROW SECURITIES
   C1v(ik2_n,ik3_n,ik4_n,  is1_n, is2_n) =C1
   L1v(ik2_n,ik3_n,ik4_n,  is1_n, is2_n) = L1
   C2v(ik2_n,ik3_n,ik4_n,  is1_n, is2_n) =C2
   L2v(ik2_n,ik3_n,ik4_n,  is1_n, is2_n) = L2
   Y1v(ik2_n,ik3_n,ik4_n,  is1_n, is2_n)= Y1
   pf_sols(ik2_n,ik3_n,ik4_n,  is1_n, is2_n)=pf

   phv(ik2_n,ik3_n,ik4_n,  is1_n, is2_n)=ph
   OmegaSplusv(ik2_n,ik3_n,ik4_n,  is1_n, is2_n)=OmegaSplus
! write(*,'(A,f15.5,2x,A,2x,f15.5,A)')ANSI_YELLOW//'Value functions at t=0: V1 = ',x0(1),'V2 = ',x0(2),ANSI_RESET


end subroutine


