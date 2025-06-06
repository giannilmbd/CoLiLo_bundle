subroutine solve_policy_function(kappa, converged)

   use globals
   use bspline_oo_module
   use interp_mod, only: bspline4_init, bspline4_eval

   use my_toolbox, only: spline_interp, fzero, toc
   use omp_lib
   use optimization_interface_mod
   use mymoments
   ! use NEQNF_INT
   ! use NEQBF_INT
   ! use NEQBJ_INT
   ! use IMSL_LIBRARIES

   implicit none
   real(wp), intent(in) :: kappa
   logical, intent(out) :: converged
   real(wp), dimension(2) :: x_both, FVEC, xguess, x_out,lb,ub
   real(wp) ::  RPARAM(5),con_lev_all(itermax)
   integer ::IPARAM(6)

   ! TYPE(gpf):: gp
   real(wp) :: parin0(100), FNORM
   logical :: decreasing
   integer :: ik2,ik3,ik4, is1, is2, trheadnum,iflag


   !include 'nlopt.f'
   ! external myfunc_FPI
   external myfunc_update

   iter = 0
   con_lev_v1 = 1e10_wp
   con_lev_v2 = 1e10_wp
   ub = [10.0_wp, 10.0_wp]
   lb = [0.009_wp, 0.009_wp]
   
   converged = .false.
   
   do while ((iter < itermax) .and. ((con_lev_v1 > sig) .or. (con_lev_v2 > sig)))
      iter = iter + 1
      ! interpolate coefficients
!$omp PARALLEL  default(shared) private(x_both,check,parin0,FVEC,lb,ub)
!$omp DO collapse(2) schedule(static)
      do is2 = 1, NS
         do is1 = 1, NS
            call spline_interp(v1(:,:,:,is1, is2), coeff_v1(:,:,:, is1, is2))
            call spline_interp(v2(:,:,:,is1, is2), coeff_v2(:,:,:, is1, is2))
         end do
      end do
!$omp END DO
      !!$omp END target teams distribute PARALLEL DO
      parin0=0.0_wp
!$omp DO   collapse(5) schedule(dynamic,1)
do is2 = 1, NS
   do is1 = 1, NS
      do ik4 = 0, NK4
         do ik3 = 0, NK3
            do ik2 = 0, NK2

               parin0(1) = real(is1, 8)
               parin0(2) = real(is2, 8)
               parin0(3) = kappa
               parin0(4) = real(ik2,8)
               parin0(5) = real(ik3,8)
               parin0(6) = real(ik4,8)

               x_both(1) = v1(ik2,ik3,ik4, is1, is2)
               x_both(2) = v2(ik2,ik3,ik4, is1, is2)

               select case (trim(inner_loop))
               case ('fzero')
                  if(iter>50)then
                     call optimize(xin=x_both, func_=solve_value, lb=lb, ub=ub, fdata=parin0, code_=12, constr=0)
                  end if 
                  call myfunc_update(2, x_both, FVEC, parin0)
               case ('update')
                  call myfunc_update(2, x_both, FVEC, parin0)
               end select

               v1_new(ik2,ik3,ik4, is1, is2) = x_both(1)
               v2_new(ik2,ik3,ik4, is1, is2) = x_both(2)

            end do
         end do
      end do
   end do
end do
!$omp END DO
!$omp END PARALLEL

      ! get convergence level
      con_lev_v1 = maxval(abs(v1_new - v1)/max(abs(v1), 1e-10_wp))
      con_lev_v2 = maxval(abs(v2_new - v2)/max(abs(v2), 1e-10_wp))
      con_lev_all(iter)=norm2([con_lev_v1, con_lev_v2],1)
      if(.false.) then
         if (iter > 5) then
            decreasing = all( con_lev_all(iter-5+1:iter)                       &
               < con_lev_all(iter-5  :iter-1) )
            if (.not. decreasing) then
               write(*,*) ANSI_YELLOW//'WARNING: Residuals are increasing. switching to cubic spline.'//ANSI_RESET
               cubic_spline = .true.
               smooth_outer = smooth_alittle
            else
               write(*,*) ANSI_GREEN//'WARNING: switching to linear interpolation.'//ANSI_RESET
               cubic_spline = .false.
               smooth_outer = smooth_alot
            end if
         end if
      end if


      !! PRINT ITERATIONS OUTCOME
      ! call toc()
      if (verbose) then
         ! write (*, '(A10,A2,i5,A2,A15,A2,f20.11,A2,A15,A2,f20.11,A2,A15,A2,f20.11,A2,A15,A2,f20.11)') 'Iteration:',',', iter,',', 'Resid ValueF 1: ',',', &
         ! con_lev_v1,',', 'Resid ValueF 2:',',', con_lev_v2,',','Mean Res v1',',',meanval_v1,',','Mean Res v2',',',meanval_v2
         write (*, '(A10,A2,i5,A2,A15,A2,e20.6,A2,A15,A2,e20.6)') 'Iteration:', ',', iter, ',', 'Resid ValueF 1: ', ',', &
            con_lev_v1, ',', 'Resid ValueF 2:', ',', con_lev_v2
         write(*,'(A,f18.8,A,f18.8,A,A,f15.7,A)')ANSI_YELLOW//'Value functions at t=0: V1 = ',v1_new(pos_zeroA1,pos_zeroA2,pos_zeroOm,pos_zero1,pos_zero2),&
            '; V2 = ',v2_new(pos_zeroA1,pos_zeroA2,pos_zeroOm,pos_zero1,pos_zero2),ANSI_RESET,ANSI_BLUE//"kappa: ",kappa,ANSI_RESET

         call toc()
      end if

      !!! check for convergence
      if ((con_lev_v1 < sig) .and. (con_lev_v2 < sig)) then
         v1 = v1_new
         v2 = v2_new 
         write(*,*) ANSI_YELLOW//'-----------------------------------------------------------'//ANSI_RESET
         write (*, '(A,i5,2x,A,es20.5,2x,A,es20.5,A)') ANSI_RED//'CONVERGED at Iteration:', &
            iter, ANSI_RESET//ANSI_BLUE//'Resid ValueF 1: ', con_lev_v1, 'Resid ValueF 2: ', con_lev_v2,ANSI_RESET
         write(*,'(A,f15.5,2x,A,2x,f15.5,A)')ANSI_YELLOW//'Value functions at t=0: V1 = ',v1(pos_zeroA1,pos_zeroA2,pos_zeroOm,pos_zero1,pos_zero2),&
            'V2 = ',v2(pos_zeroA1,pos_zeroA2,pos_zeroOm,pos_zero1,pos_zero2),ANSI_RESET
         call toc()
         write(*,*) ANSI_YELLOW//'-----------------------------------------------------------'//ANSI_RESET
         
         converged = .true.
      else
         v1 = smooth_outer*v1_new + (1e0_wp - smooth_outer)*v1
         v2 = smooth_outer*v2_new + (1e0_wp - smooth_outer)*v2
        ! the interpolation is at the top of the loop    
      end if !if converged

   end do ! iterations
!$omp PARALLEL  default(shared)  ! must interpolate after solution (too; or only in this special case)
!$omp DO collapse(2) schedule(static)
   do is2 = 1, NS
      do is1 = 1, NS

         call spline_interp(v1_new(:,:,:, is1, is2), coeff_v1(:,:,:, is1, is2))
         call spline_interp(v2_new(:,:,:, is1, is2), coeff_v2(:,:,:, is1, is2))


      end do
   end do
   !!$omp END target teams distribute PARALLEL DO
   !$omp END PARALLEL

end subroutine solve_policy_function
