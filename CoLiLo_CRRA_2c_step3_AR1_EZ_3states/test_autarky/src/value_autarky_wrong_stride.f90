subroutine value_autarky()

   use globals
   use my_toolbox
   use omp_lib

   implicit none
   real(kind=wp) :: v1plus, v2plus, C1, C2, Ev1aut, Ev2aut, C1plus,ph,Q,x_states(2)
   real(kind=wp) :: A1plus1, A2plus1, A1_n, A2_n, pf, pf_new, L1, L2,  Y1, res, pred, parin01(4)
   logical :: check_
   integer ::  ik2,ik3,is2, is_n2, is_p2, is1, is_n1, is_p1, cnt
   ! start from solution of CM

   iter = 0
   con_lev_v1 = 1e10_wp
   con_lev_v2 = 1e10_wp
   parin01 = 0.0_wp
   outer_loop: do while ((iter < itermax) .and. ((con_lev_v1 > sig) .or. (con_lev_v2 > sig)))
      iter = iter + 1
      do is_n1 = 1, NS
        do is_n2 = 1, NS

           call spline_interp(v1_aut(:, :, is_n1, is_n2), coeff_v1_aut(:, :, is_n1, is_n2))
           call spline_interp(v2_aut(:, :, is_n1, is_n2), coeff_v2_aut(:, :, is_n1, is_n2))
        end do
     end do


!$omp PARALLEL default(shared)
      !$omp DO private(Ev1aut,Ev2aut,pf,C1,L1,C2,L2,A1plus1,A2plus1,ph,Q,x_states,v1plus,v2plus) &
      !$omp     collapse(4) schedule(dynamic,1)
      do ik2 = 0,NK2
         do ik3 = 0,NK3
            do is_n1 = 1, NS
               do is_n2 = 1, NS
                  pf=minval([maxval([price_down,pf_sols_aut(ik2,ik3,is_n1,is_n2)]),price_up],1);
                  A1plus1 = A1(ik2)**rho_A1* exp(eta1(is_n1))
                  A2plus1 = A2(ik3)**rho_A2* exp(eta2(is_n2))
                  call bisection_price_autarky(pf, C1, L1, C2, L2, A1plus1, A2plus1,ph,Q)
                  pf_sols_aut(ik2,ik3,is_n1,is_n2)=pf
                  Ev1aut = 0e0_wp
                  Ev2aut = 0e0_wp
                  
                  
                  x_states=[min(A1plus1, UB_all(2)), min(A2plus1, UB_all(3))]
                  do is_p1 = 1, NS
                     do is_p2 = 1, NS
                        v1plus = (spline_eval(x_states, &
                           coeff_v1_aut(:, :, is_p1, is_p2), LB_all(2:3), UB_all(2:3)))
                        v2plus = (spline_eval(x_states, &
                           coeff_v2_aut(:, :, is_p1, is_p2), LB_all(2:3), UB_all(2:3)))
                        Ev1aut = Ev1aut + pi1(is_n1, is_p1)*pi2(is_n2, is_p2)*(v1plus**theta) ! See Rudebush and Swansson on changing sign with 1/gamma>1
                        Ev2aut = Ev2aut + pi1(is_n1, is_p1)*pi2(is_n2, is_p2)*(v2plus**theta)

                     end do
                  end do

                  v1_aut_new(ik2,ik3,is_n1,is_n2)  = -(1e0_wp - beta)*(C1**(1e0_wp - 1e0_wp/gamma)/(1e0_wp - 1e0_wp/gamma)-chhi*L1**(1.0e0_wp+varphi)/(1.0e0_wp+varphi)) +&
                   beta*Ev1aut**(1.0_wp/theta) ! check sidn Rud&Swans
                  v2_aut_new(ik2,ik3,is_n1,is_n2)  = -(1e0_wp - beta)*(C2**(1e0_wp - 1e0_wp/gamma)/(1e0_wp - 1e0_wp/gamma)-chhi*L2**(1.0e0_wp+varphi)/(1.0e0_wp+varphi)) +&
                   beta*Ev2aut**(1.0_wp/theta)
                  C1v_aut(ik2,ik3,is_n1,is_n2) = C1
                  C2v_aut(ik2,ik3,is_n1,is_n2) = C2
                  L1v_aut(ik2,ik3,is_n1,is_n2) = L1
                  L2v_aut(ik2,ik3,is_n1,is_n2) = L2
                  Y1v_aut(ik2,ik3,is_n1,is_n2) = A1plus1*L1**(1.0_wp-alphha)
                  Y2v_aut(ik2,ik3,is_n1,is_n2) = A2plus1*L2**(1.0_wp-alphha)
                  Qv_aut(ik2,ik3,is_n1,is_n2) = Q
                  pfv_aut(ik2,ik3,is_n1,is_n2) = pf
                  phv_aut(ik2,ik3,is_n1,is_n2) = ph
               end do
            end do
         end do
      end do

!$omp END PARALLEL
   ! get convergence level
   con_lev_v1 = maxval(abs(v1_aut_new - v1_aut)/max(abs(v1_aut), 1e-10_wp))
   con_lev_v2 = maxval(abs(v2_aut_new - v2_aut)/max(abs(v2_aut), 1e-10_wp))
! con_lev_s1 = maxval(abs(s1_new - s1)/max(abs(s1), 1d-10))

!! PRINT ITERATIONS OUTCOME
! call toc()
   if (verbose) then
      write (*, '(A30,i5,2x,A15,e20.11,2x,A15,e20.11)') ANSI_YELLOW//'Iteration AUTARKY: '//ANSI_RESET, iter, 'Resid ValueF 1: ', &
         con_lev_v1, 'Resid ValueF 2:', con_lev_v2
   end if

!!! check for convergence
   if ((con_lev_v1 < sig) .and. (con_lev_v2 < sig)) then
      !call output(kappa) ! THIS CHECKS GRAPHICALLY THE RESULTS AND COMPUTES EULER RESIDUAL
      write(*,'(A,i5,2x,A15,es20.5,2x,A15,es20.5)') '\033[95m AUTARKY CONVERGED at Iteration: \033[0m', iter, '[34m Resid ValueF 1: [0m',con_lev_v1, '[34m Resid ValueF 2: [0m', con_lev_v2

      ! return
   else

      v1_aut = smooth_aut*v1_aut_new + (1e0_wp - smooth_aut)*v1_aut
      v2_aut = smooth_aut*v2_aut_new + (1e0_wp - smooth_aut)*v2_aut
   end if

end do outer_loop ! iterations
end subroutine
