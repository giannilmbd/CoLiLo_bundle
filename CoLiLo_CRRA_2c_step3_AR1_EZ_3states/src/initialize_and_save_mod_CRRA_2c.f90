module initialize_and_save_mod_CRRA_2c
   use my_kinds_mod, only: wp
   use globals
   implicit none
contains

   subroutine initialvalues_CRRA(dir,cntry, def_cntry)
      use omp_lib
      use my_I_O
      use globals
      use read_csv_mod
      use my_kinds_mod
      use optimization_interface_mod
      implicit none

      integer ::  ik2,ik3,ik4, is1, is2
      real(kind=wp) :: Omega1, Omega2, C1, C2,mu, OmegaSplus, ph, pf, pC, L1, L2, Q, Y1, parin0(100)
      real(kind=wp) ::res, pred, kappa, lb(1), ub(1), pfdummy(1),A1plus1,A2plus1
      character(len=*), intent(in) :: cntry, def_cntry
      character(len=100) :: def_filename, filename,dir
      logical :: file_exists, def_file_exists
      pf_sols = 1.0e0_wp

      if (doinitval) then
         !$omp PARALLEL DO private(Omega1,Omega2,OmegaSplus,C1, C2, L1, L2, Y1, Q, res, pred, ph, pf, pfdummy, &
         !$omp&                       kappa, lb, ub, parin0, A1plus1, A2plus1) &
         !$omp&     default(shared) &
         !$omp&     collapse(5)    &
         !$omp&     schedule(dynamic,1)

         
            do ik2 =0, NK2
               do ik3 =0, NK3
                  do ik4 = 0, NK4
                  do is1 = 1, NS
                     do is2 = 1, NS

                        pfdummy(1) = pf_sols(ik2,ik3,ik4, is1, is2)
                        kappa = 1.0_wp 
                        lb = price_down
                        ub = price_up
                        parin0(1) = real(is1, 8)
                        parin0(2) = real(is2, 8)
                        parin0(3) = kappa
                        parin0(4) = real(ik2, 8)
                        parin0(5) = real(ik3, 8)
                        parin0(6) = real(ik4, 8)
                        A1plus1= A1(ik2)**rho_A1*exp(eta1(is1))
                        A2plus1= A2(ik3)**rho_A2*exp(eta2(is2))
                        Omega1 = beta*(A1(ik2))**(1e0_wp/theta - 1e0_wp)*(A1plus1**(1e0_wp - 1e0_wp/gamma)/(1e0_wp - 1e0_wp/gamma))**(theta - 1)

                        Omega2 = beta*(A2(ik3))**(1e0_wp/theta - 1e0_wp)*(A2plus1**(1e0_wp - 1e0_wp/gamma)/(1e0_wp - 1e0_wp/gamma))**(theta - 1)
                                           ! update Omega
                                           OmegaSplus =Omega1/Omega2* OmegaS(ik4)

                        !   ! Start the timer
                        ! call CPU_TIME(start_time)
                        ! call optimize(xin=pfdummy, func_=solve_pf, lb=lb, ub=ub, fdata=parin0, code_=12)
                        ! pf = pfdummy(1)
                        call core_model(pf=1.0_wp, kappa=kappa,  OmegaS_in=OmegaSplus, C=C1, Cs=C2, L=L1, Ls=L2, &
                           ph=ph, Q=Q, Y1=Y1, res=res, pred=pred, D=A1plus1, Ds=A2plus1)
                        ! future X state
                        pf_sols(ik2,ik3,ik4, is1, is2) = pf

                        ! note that preferences are V=(1-beta)*U(c,l)+beta*E(V)
                        v1(ik2,ik3,ik4, is1, is2) = -(C1**(1e0_wp - 1e0_wp/gamma)/(1e0_wp - 1e0_wp/gamma) - &
                           chhi*L1**(1.0e0_wp + varphi)/(1.0e0_wp + varphi))  ! check sidn Rud&Swans
                        v2(ik2,ik3,ik4, is1, is2) = -(C2**(1e0_wp - 1e0_wp/gamma)/(1e0_wp - 1e0_wp/gamma) - &
                           chhi*L2**(1.0e0_wp + varphi)/(1.0e0_wp + varphi))

                        s1(ik2,ik3,ik4, is1, is2) = (C1**(-1e0_wp/gamma)*(C1 - ph*Y1))/(1.0_wp - beta)

                     end do
                  end do
               end do
            end do
         end do
      
         !$omp END PARALLEL DO

         pf_sols_aut = pf_sols(:,:,pos_zeroOm, :, :)
         v1_aut = v1(:,:,pos_zeroOm, :, :);
         v2_aut = v2(:,:,pos_zeroOm, :, :);

      else ! read from saved policy functions

         ! fist check if policy functions exist either for each country or for default country
         ! defaulting to def_cntry if no policy function was stored for specific country

         filename = trim(StrStrStr(dir=trim(dir), str1='v1_data_'//mom2loop//'_', cnt=cntry, str2='.csv'))
         inquire (file=filename, exist=file_exists)
         def_filename = trim(StrStrStr(dir=trim(dir), str1='v1_data_'//mom2loop//'_', cnt=def_cntry, str2='.csv'))
         inquire (file=def_filename, exist=def_file_exists)
         if (file_exists) then
            write (*, *) ANSI_BLUE//'FOUND FILES FOR POLICY FUCNTION COUNTRY COUNTRY ==>>'//cntry//ANSI_RESET
            call read_csv(filename=trim(filename), data=matrix_csv_v1)
            call read_csv(filename=trim(StrStrStr(dir=trim(dir), str1='v2_data_'//mom2loop//'_', cnt=cntry, str2='.csv')), data=matrix_csv_v2)

            call read_csv(filename=trim(StrStrStr(dir=trim(dir), str1='s1_data_'//mom2loop//'_', cnt=cntry, str2='.csv')), data=matrix_csv_s1)

            call read_csv(filename=trim(StrStrStr(dir=trim(dir), str1='pB_data_'//mom2loop//'_', cnt=cntry, str2='.csv')), data=matrix_csv_pf)

            ! # AUTARKY
            call read_csv(filename=trim(StrStrStr(dir=trim(dir),str1='v1_aut_data_'//mom2loop//'_', cnt=cntry, str2='.csv')), data=matrix_csv_aut_v1)
            call read_csv(filename=trim(StrStrStr(dir=trim(dir),str1='v2_aut_data_'//mom2loop//'_', cnt=cntry, str2='.csv')), data=matrix_csv_aut_v2)

            call read_csv(filename=trim(StrStrStr(dir=trim(dir), str1='pB_aut_data_'//mom2loop//'_', cnt=cntry, str2='.csv')), data=matrix_csv_aut_pf)

         elseif (def_file_exists) then ! USE A DEFAULT 'COUNTRY' (first of the list); if even this does not exist will stop with an error
            write (*, *) ANSI_RED//'NOT FOUND FILES FOR POLICY COUNTRY ==>'//cntry
            write (*,*)  ANSI_RED//'LOADING DEFAULT COUNTRY ==>>'//def_cntry//ANSI_RESET

            call read_csv(filename=trim(def_filename), data=matrix_csv_v1)
            call read_csv(filename=trim(StrStrStr(dir=trim(dir), str1='v2_data_'//mom2loop//'_', cnt=def_cntry, str2='.csv')), data=matrix_csv_v2)

            call read_csv(filename=trim(StrStrStr(dir=trim(dir), str1='s1_data_'//mom2loop//'_', cnt=def_cntry, str2='.csv')), data=matrix_csv_s1)

            call read_csv(filename=trim(StrStrStr(dir=trim(dir), str1='pB_data_'//mom2loop//'_', cnt=def_cntry, str2='.csv')), data=matrix_csv_pf)

            ! # AUTARKY
            call read_csv(filename=trim(StrStrStr(dir=trim(dir),str1='v1_aut_data_'//mom2loop//'_', cnt=def_cntry, str2='.csv')), data=matrix_csv_aut_v1)
            call read_csv(filename=trim(StrStrStr(dir=trim(dir),str1='v2_aut_data_'//mom2loop//'_', cnt=def_cntry, str2='.csv')), data=matrix_csv_aut_v2)

            call read_csv(filename=trim(StrStrStr(dir=trim(dir), str1='pB_aut_data_'//mom2loop//'_', cnt=def_cntry, str2='.csv')), data=matrix_csv_aut_pf)

         else
            write (*, *) ANSI_RED//' NO INITIAL CONDITION FILES WERE FOUND EVEN FOR DEFAULT COUNTRY ==> '//def_cntry//ANSI_RESET
            stop
         end if

         v1 = RESHAPE(source=(matrix_csv_v1), shape=shape(v1))
         v2 = RESHAPE(source=(matrix_csv_v2), shape=shape(v2))
         s1 = RESHAPE(source=(matrix_csv_s1), shape=shape(s1))

         pF_sols = RESHAPE(source=(matrix_csv_pf), shape=shape(pf_sols))

         v1_aut = RESHAPE(source=(matrix_csv_aut_v1), shape=shape(v1_aut))
         v2_aut = RESHAPE(source=(matrix_csv_aut_v2), shape=shape(v2_aut))

         pf_sols_aut = RESHAPE(source=(matrix_csv_aut_pf), shape=shape(pf_sols_aut))

      end if !doinitval

   end subroutine initialvalues_CRRA

   subroutine save_results_CRRA(dir,cnt_mom)
      use globals,dummy=>write_csv_portable
      use my_I_O
      use read_csv_mod
      use write_csv_portable_module,only: write_csv_portable
      implicit none

      character(len=*), intent(in) :: cnt_mom
      character(len=100) :: filename, dir2, direnq,dir
      integer ::  unitnumber1, unitnumber2, unitnumber3, nlines,status
      logical :: dir_exists
      ! 101 format(1x, *(g0, ", "))
      ! SAVE THE SOLUTION TO FILE values and coefficients

      dir2 = trim(dir)
      ! inquire (directory=trim(dir2), exist=dir_exists)
      call execute_command_line('test -d "' // trim(dir2) // '"', exitstat=status)
      dir_exists = (status == 0)
      if (dir_exists .neqv. .true.) then
         direnq = 'mkdir -p'//trim(dir2)
         call execute_command_line(direnq)
      end if

      matrix_csv_v1 = RESHAPE(source=v1, &
         shape=shape(matrix_csv_v1))
      matrix_csv_v2 = RESHAPE(source=v2, &
         shape=shape(matrix_csv_v2))
      matrix_csv_s1 = RESHAPE(source=s1, &
         shape=shape(matrix_csv_s1))

      matrix_csv_pf = RESHAPE(source=pf_sols, &
         shape=shape(matrix_csv_pf))
      filename = trim(StrStrStr(dir=trim(dir), str1='v1_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv'))
      call write_csv_portable(filename=filename, data=matrix_csv_v1)
      filename = trim(StrStrStr(dir=trim(dir), str1='v2_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv'))
      call write_csv_portable(filename=filename, data=matrix_csv_v2)
      filename = trim(StrStrStr(dir=trim(dir), str1='s1_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv'))
      call write_csv_portable(filename=filename, data=matrix_csv_s1)
      filename = trim(StrStrStr(dir=trim(dir), str1='pB_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv'))
      call write_csv_portable(filename=filename, data=matrix_csv_pf)

      ! SAVE AUTARKY
      matrix_csv_aut_v1 = RESHAPE(source=v1_aut, &
         shape=shape(matrix_csv_aut_v1))
      matrix_csv_aut_v2 = RESHAPE(source=v2_aut, &
         shape=shape(matrix_csv_aut_v2))

      matrix_csv_aut_pf = RESHAPE(source=pf_sols_aut, &
         shape=shape(matrix_csv_aut_pf))
      filename = trim(StrStrStr(dir=trim(dir), str1='v1_aut_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv'))
      call write_csv_portable(filename=filename, data=matrix_csv_aut_v1)
      filename = trim(StrStrStr(dir=trim(dir), str1='v2_aut_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv'))
      call write_csv_portable(filename=filename, data=matrix_csv_aut_v2)
      filename = trim(StrStrStr(dir=trim(dir), str1='pB_aut_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv'))
      call write_csv_portable(filename=filename, data=matrix_csv_aut_pf)
      ! save coefficients
      matrix_coefs_csv_v1 = RESHAPE(source=coeff_v1, &
         shape=shape(matrix_coefs_csv_v1))
      matrix_coefs_csv_v2 = RESHAPE(source=coeff_v2, &
         shape=shape(matrix_coefs_csv_v2))

      matrix_coefs_csv_s1 = RESHAPE(source=coeff_s1, &
         shape=shape(matrix_coefs_csv_s1))
      filename = trim(StrStrStr(dir=trim(dir), str1='coef_v1_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv'))
      call write_csv_portable(filename=filename, data=matrix_coefs_csv_v1)
      filename = trim(StrStrStr(dir=trim(dir), str1='coef_v2_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv'))
      call write_csv_portable(filename=filename, data=matrix_coefs_csv_v2)
      filename = trim(StrStrStr(dir=trim(dir), str1='coef_s1_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv'))
      call write_csv_portable(filename=filename, data=matrix_coefs_csv_s1)

      ! TEST IF SAVING LOSES DATA
      if (.false.) then !mute if not debugging
         nlines = size(matrix_csv_v1, 1)
         filename = StrStrStr(dir=trim(dir), str1='v1_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv')
         call read_csv(filename=filename, data=matrix_csv_v1, nlines=nlines)

         v1_test = RESHAPE(source=(matrix_csv_v1), shape=shape(v1_test))
         meanval_v1 = maxval(abs(v1_test - v1))
         print *, '\033[34m MAX DISCREPANCY RELATIVE TO SAVED DATA \033[0m', meanval_v1

      end if

   end subroutine

end module
