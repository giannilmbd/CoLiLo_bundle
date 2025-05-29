subroutine solve_arrow_securities()
   use my_kinds_mod, only: wp,sp
   use my_toolbox, only: spline_eval,spline_interp,toc
   use globals

   implicit none
   ! Declare the arguments
  
   integer :: ik2, ik3, ik4, is1, is2,is_p1,is_p2
   real(kind=wp) :: Omega1plus,s1plus,x_states(3)
   real(kind=wp) :: A1plus1_precomputed(0:NK2,NS), A2plus1_precomputed(0:NK3,NS),A1plus1,A2plus1

   iter=0
   con_lev_s1 = 10.0e6_wp

   ! Precompute A1plus1 values using is1 (shock to country 1)
   do is1 = 1, NS
      do ik2 = 0, NK2
         A1plus1_precomputed(ik2, is1) = A1(ik2)**rho_A1 * exp(eta1(is1))
      end do
   end do

   ! Precompute A2plus1 values using is2 (shock to country 2)
   do is2 = 1, NS
      do ik3 = 0, NK3
         A2plus1_precomputed(ik3, is2) = A2(ik3)**rho_A2 * exp(eta2(is2))
      end do
   end do

   iteration: do while ((iter < itermax) .and. ((con_lev_s1 > sig_arrow)))
      iter = iter + 1
!$omp PARALLEL  default(shared)
      !$omp DO collapse(2) schedule(static)
      do is2 = 1, NS
         do is1 = 1, NS
            call spline_interp(s1(:,:,:,is1, is2), coeff_s1(:,:,:,  is1, is2))
         end do
      end do
      !$omp END DO
      !$omp END PARALLEL

      Es1=0.0_wp
      ! Loop over the grid points for the Arrow securities
      !$omp PARALLEL default(shared) private(Omega1plus,s1plus,A1plus1,A2plus1,x_states)
      !$omp DO collapse(5) schedule(static)
      do is2 = 1, NS
         do is1 = 1, NS
            do ik4 = 0, NK4
               do ik3 = 0, NK3
                  do ik2 = 0, NK2
                     ! form expected value of Arrow security
                     A1plus1 = A1plus1_precomputed(ik2, is1)  ! Using is1 for country 1
                     A2plus1 = A2plus1_precomputed(ik3, is2)  ! Using is2 for country 2

                     x_states=[A1plus1, A2plus1,OmegaSplusv(ik2,ik3,ik4,is1,is2)]
                     do is_p2 = 1, NS
                        do is_p1 = 1, NS
                           s1plus = spline_eval(x_states,coeff_s1( :, :,:, is_p1, is_p2), LB_all, UB_all)
                           Es1(ik2, ik3,ik4,is1, is2) = Es1(ik2, ik3,ik4,is1, is2) + &
                              pi1(is1, is_p1)*pi2(is2, is_p2)*(v1(ik2,ik3,ik4,is_p1,is_p2))**(theta - 1)*s1plus
                        end do
                     end do
                     s1_new(ik2,ik3,ik4,is1,is2) = C1v(ik2,ik3,ik4,is1,is2)**(-1.0_wp/gamma)*&
                        (C1v(ik2,ik3,ik4,is1,is2)-phv(ik2,ik3,ik4,is1,is2)*Y1v(ik2,ik3,ik4,is1,is2))+&
                        beta*(Ev1(ik2, ik3,ik4,is1,is2))**(1.0_wp/theta - 1.0_wp)*Es1(ik2,ik3,ik4,is1, is2)

                  end do
               end do
            end do
         end do
      end do
!$omp END DO
      !$omp END PARALLEL
      con_lev_s1 = maxval(abs(s1_new - s1))

      if ((con_lev_s1 < sig_arrow)) then
         write (*, '(A,i5,2x,A15,es20.5,A)') ANSI_BLUE//'CONVERGED at Iteration: ', &
            iter, 'Resid Arrow Securities:', con_lev_s1,ANSI_RESET
         write(*,'(A,f15.5,2x,A,f15.5,A)')  ANSI_GREEN//'Arrow 1: '//ANSI_RED//'MAX = ', maxval((s1_new)),', MIN = ',minval((s1_new)),ANSI_RESET

         call toc()
         s1 = s1_new
      else
         if(verbose) then
            write (*, '(A,i5,2x,A15,es20.5,2x,A,f15.8,A)') ANSI_BLUE//'Iteration: ', &
               iter, 'Resid Arrow Securities:', con_lev_s1,ANSI_RESET//ANSI_RED//"Arrow at t=0: ",s1(pos_zeroA1,pos_zeroA2,pos_zeroOm,pos_zero1,pos_zero2),ANSI_RESET
         end if
         s1 = smooth_arrow*s1_new + (1e0_wp - smooth_arrow)*s1
      end if

   end do iteration

   ! !$omp PARALLEL  default(shared)
   !    !$omp DO collapse(2) schedule(static)
   ! do is2 = 1, NS
   !    do is1 = 1, NS
   !       call spline_interp(s1(:,:,:,is1, is2), coeff_s1(:,:,:,  is1, is2))
   !    end do
   ! end do
   ! !$omp END DO
   ! !$omp END PARALLEL

end subroutine solve_arrow_securities

