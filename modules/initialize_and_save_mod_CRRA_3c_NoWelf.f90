module initialize_and_save_mod_CRRA_3c_NoWelf
    use my_kinds_mod, only: wp
    use globals

    implicit none
contains

    subroutine initialvalues_CRRA_3c(dir,cntry,def_cntry)
        use omp_lib
        use my_I_O
        use globals
        use read_csv_mod
        use my_kinds_mod
        use optimization_interface_mod
        use minpack_module, only: hybrd1, dpmpar, enorm

        implicit none

        integer :: ik1, ik2, is1, is2, is3,info
        integer , parameter :: n_prices=2
        real(kind=wp) :: C1, C2, C3, Yworld, pA, pB, pC, L1, L2, L3, Qb, Qc, Y1, Y2, parin0(7)
        real(kind=wp), dimension(n_prices) ::res, pred, pf, kappa, lb, ub
        character(len=*), intent(in) :: cntry,def_cntry,dir

        character(len=100) :: def_filename,filename
        logical :: file_exists,def_file_exists
        integer, parameter :: lwa = (n_prices*(3*n_prices + 13))/2
        real(kind=wp) :: wa(lwa)

! pB_sols=1.0e0_wp
! pC_sols=1.0e0_wp
        if (doinitval) then
            !$omp PARALLEL DO  private(C1,C2,C3,L1,L2,L3,Y1,Y2,Qb,Qc,res,pred,pA,pf,kappa,lb,ub,parin0) default(shared) collapse(5)

            do ik1 = 0, NK1
                do ik2 = 0, NK2
                    do is1 = 1, NS
                        do is2 = 1, NS
                            do is3 = 1, NS
                                pf(1)=pB_sols(ik1,ik2,is1,is2,is3)
                                pf(2)=pC_sols(ik1,ik2,is1,is2,is3)
                                kappa(1) = kappaB(ik1)
                                kappa(2) = kappaC(ik2)
                                lb = price_down
                                ub = price_up
                                parin0(1) = real(is1, wp)
                                parin0(2) = real(is2, wp)
                                parin0(3) = real(is3, wp)
                                parin0(4) = real(ik1, wp)
                                parin0(5) = real(ik2, wp)
                                !   ! Start the timer
                                ! call CPU_TIME(start_time)
                                ! call optimize(xin=pf, func_=solve_pf, lb=lb, ub=ub, fdata=parin0, code_=25)
                                call hybrd1(solve_pf_minpack, size(pf,1), pf, res, sig_pf, info, wa, lwa, parin0)
                                call core_model(pf=pf, kappa=kappa, C=C1, Cs=C2, Cc=C3, L=L1, Ls=L2, Lc=L3, pA=pA, &
                             Qb=Qb, Qc=Qc, Ya=Y1, Yb=Y2, res=res, pred=pred, D=exp(eta1(is1)), Ds=exp(eta2(is2)), Dc=exp(eta3(is3)))
                                ! future X state
                             pB_sols(ik1,ik2,is1,is2,is3)=pf(1)
                             pC_sols(ik1,ik2,is1,is2,is3)=pf(2)
                             ! note that preferences are V=(1-beta)*U(c,l)+beta*E(V)
                                v1(ik1, ik2, is1, is2, is3) = (C1**(1e0_wp - 1e0_wp/gamma)/(1e0_wp - 1e0_wp/gamma) - &
                                                              chhi*L1**(1.0e0_wp + varphi)/(1.0e0_wp + varphi))  ! check sidn Rud&Swans
                                ! Note that Arrow securities are such that A*C**(-1/gamma)=C**(-1/gamma)*(C-p*Y)+beta*E(C**(-1/gamma)A)
                                s1(ik1, ik2, is1, is2, is3) = (C1**(-1e0_wp/gamma)*(C1 - pA*Y1))/(1.0_wp-beta)
                                s2(ik1, ik2, is1, is2, is3) = (C2**(-1e0_wp/gamma)*(C2 - pf(1)/Qb*Y2))/(1.0_wp-beta)
                            end do
                        end do
                    end do

                end do
            end do
            !$omp END PARALLEL DO
            pB_sols_aut=pB_sols(pos_zeroKB, pos_zeroKC, :, :, :);
            pC_sols_aut=pC_sols(pos_zeroKB, pos_zeroKC, :, :, :);
            v1_aut = v1(pos_zeroKB, pos_zeroKC, :, :, :); 
        else

 ! fist check if policy functions exist either for each country or for default country            
            ! defaulting to def_cntry if no policy function was stored for specific country

            filename=trim(dir)//'/v1_data_'//trim(mom2loop)//'_'//trim(cntry)//'.csv'
            inquire(file=trim(filename), exist=file_exists)          
            def_filename=trim(StrStrStr(dir=trim(dir), str1='v1_data_'//mom2loop//'_', cnt=def_cntry, str2='.csv'))
      inquire(file=def_filename, exist=def_file_exists)
      if (file_exists) then
      write(*,*) ANSI_BLUE//'FOUND FILES FOR POLICY FUCNTION COUNTRY ==>>'//cntry//ANSI_RESET
      call read_csv(filename=trim(filename), data=matrix_csv_v1)
      call read_csv(filename=trim(StrStrStr(dir=trim(dir), str1='s1_data_'//mom2loop//'_', cnt=cntry, str2='.csv')), data=matrix_csv_s1)
      call read_csv(filename=trim(StrStrStr(dir=trim(dir), str1='s2_data_'//mom2loop//'_', cnt=cntry, str2='.csv')), data=matrix_csv_s2)

      call read_csv(filename=trim(StrStrStr(dir=trim(dir), str1='pB_data_'//mom2loop//'_', cnt=cntry, str2='.csv')), data=matrix_csv_pB)
      call read_csv(filename=trim(StrStrStr(dir=trim(dir), str1='pC_data_'//mom2loop//'_', cnt=cntry, str2='.csv')), data=matrix_csv_pC)

                  ! # AUTARKY
      call read_csv(filename=trim(StrStrStr(dir=trim(dir),str1='v1_aut_data_'//mom2loop//'_', cnt=cntry, str2='.csv')), data=matrix_csv_aut_v1)
      
      call read_csv(filename=trim(StrStrStr(dir=trim(dir), str1='pB_aut_data_'//mom2loop//'_', cnt=cntry, str2='.csv')), data=matrix_csv_aut_pB)
      call read_csv(filename=trim(StrStrStr(dir=trim(dir), str1='pC_aut_data_'//mom2loop//'_', cnt=cntry, str2='.csv')), data=matrix_csv_aut_pC)
      
      
      
      elseif (def_file_exists) then ! USE A DEFAULT 'COUNTRY' (first of the list); if even this does not exist will stop with an error
      write(*,*) ANSI_RED//'NOT FOUND FILES FOR POLICY COUNTRY ==>'//cntry//': LOADING DEFAULT COUNTRY ==>>'//def_cntry//ANSI_RESET
      
        call read_csv(filename=trim(def_filename), data=matrix_csv_v1)
        call read_csv(filename=trim(StrStrStr(dir=trim(dir), str1='s1_data_'//mom2loop//'_', cnt=def_cntry, str2='.csv')), data=matrix_csv_s1)
        call read_csv(filename=trim(StrStrStr(dir=trim(dir), str1='s2_data_'//mom2loop//'_', cnt=def_cntry, str2='.csv')), data=matrix_csv_s2)
  
        call read_csv(filename=trim(StrStrStr(dir=trim(dir), str1='pB_data_'//mom2loop//'_', cnt=def_cntry, str2='.csv')), data=matrix_csv_pB)
        call read_csv(filename=trim(StrStrStr(dir=trim(dir), str1='pC_data_'//mom2loop//'_', cnt=def_cntry, str2='.csv')), data=matrix_csv_pC)
  
                    ! # AUTARKY
        call read_csv(filename=trim(StrStrStr(dir=trim(dir),str1='v1_aut_data_'//mom2loop//'_', cnt=def_cntry, str2='.csv')), data=matrix_csv_aut_v1)
  
        call read_csv(filename=trim(StrStrStr(dir=trim(dir), str1='pB_aut_data_'//mom2loop//'_', cnt=def_cntry, str2='.csv')), data=matrix_csv_aut_pB)
        call read_csv(filename=trim(StrStrStr(dir=trim(dir), str1='pC_aut_data_'//mom2loop//'_', cnt=def_cntry, str2='.csv')), data=matrix_csv_aut_pC)
      else
        write(*,*) ANSI_RED//' NO INITIAL CONDITION FILES WERE FOUND EVEN FOR DEFAULT COUNTRY ==> '//def_cntry//ANSI_RESET
        stop
      endif 

            v1 = RESHAPE(source=(matrix_csv_v1), shape=shape(v1))

            s1 = RESHAPE(source=(matrix_csv_s1), shape=shape(s1))
            s2 = RESHAPE(source=(matrix_csv_s2), shape=shape(s2))
            pB_sols = RESHAPE(source=(matrix_csv_pB), shape=shape(pB_sols))
            pC_sols = RESHAPE(source=(matrix_csv_pC), shape=shape(pC_sols))



            v1_aut = RESHAPE(source=(matrix_csv_aut_v1), shape=shape(v1_aut))


            pB_sols_aut = RESHAPE(source=(matrix_csv_aut_pB), shape=shape(pB_sols_aut))
            pC_sols_aut = RESHAPE(source=(matrix_csv_aut_pC), shape=shape(pC_sols_aut))
        end if !doinitval

    end subroutine initialvalues_CRRA_3c

    
subroutine save_results_CRRA_3c(dir,cnt_mom)
        use globals
        use my_I_O
        use read_csv_mod
        implicit none

        character(len=*), intent(in) :: cnt_mom
        character(len=100):: filename, dir2, direnq,dir
        integer :: ik1, unitnumber1, unitnumber2, unitnumber3, nlines
         logical :: dir_exists


          dir2 = trim(dir)
                inquire (directory=trim(dir2), exist=dir_exists)
                if (dir_exists .neqv. .true.) then
                    direnq = 'mkdir '//trim(dir2)
                    call execute_command_line(direnq)
                end if

        ! 101 format(1x, *(g0, ", "))
        ! SAVE THE SOLUTION TO FILE values and coefficients
        matrix_csv_v1 = RESHAPE(source=v1, &
                                shape=shape(matrix_csv_v1))

        matrix_csv_s1 = RESHAPE(source=s1, &
                                shape=shape(matrix_csv_s1))
        matrix_csv_s2 = RESHAPE(source=s2, &
                                shape=shape(matrix_csv_s2))

        matrix_csv_pB = RESHAPE(source=pB_sols, &
                                shape=shape(matrix_csv_pB))
        matrix_csv_pC = RESHAPE(source=pC_sols, &
                                shape=shape(matrix_csv_pC))

   call write_csv(filename=trim(StrStrStr(dir=trim(dir), str1='v1_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv')), data=matrix_csv_v1)
   call write_csv(filename=trim(StrStrStr(dir=trim(dir), str1='s1_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv')), data=matrix_csv_s1)
   call write_csv(filename=trim(StrStrStr(dir=trim(dir), str1='s2_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv')), data=matrix_csv_s2)
   call write_csv(filename=trim(StrStrStr(dir=trim(dir), str1='pB_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv')), data=matrix_csv_pB)
   call write_csv(filename=trim(StrStrStr(dir=trim(dir), str1='pC_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv')), data=matrix_csv_pC)

        ! SAVE AUTARKY
        matrix_csv_aut_v1 = RESHAPE(source=v1_aut, &
                                    shape=shape(matrix_csv_aut_v1))


        matrix_csv_aut_pB = RESHAPE(source=pB_sols_aut, &
                                shape=shape(matrix_csv_aut_pB))
        matrix_csv_aut_pC = RESHAPE(source=pC_sols_aut, &
                                shape=shape(matrix_csv_aut_pC))

        call write_csv(filename=trim(StrStrStr(dir=trim(dir),str1='v1_aut_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv')), data=matrix_csv_aut_v1)
        call write_csv(filename=trim(StrStrStr(dir=trim(dir), str1='pB_aut_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv')), data=matrix_csv_aut_pB)
        call write_csv(filename=trim(StrStrStr(dir=trim(dir), str1='pC_aut_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv')), data=matrix_csv_aut_pC)
        ! save coefficients
        matrix_coefs_csv_v1 = RESHAPE(source=coeff_v1, &
                                      shape=shape(matrix_coefs_csv_v1))


        matrix_coefs_csv_s1 = RESHAPE(source=coeff_s1, &
                                      shape=shape(matrix_coefs_csv_s1))
        matrix_coefs_csv_s2 = RESHAPE(source=coeff_s2, &
                                      shape=shape(matrix_coefs_csv_s2))

        call write_csv(filename=trim(StrStrStr(dir=trim(dir),str1='coef_v1_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv')), data=matrix_coefs_csv_v1)
        call write_csv(filename=trim(StrStrStr(dir=trim(dir),str1='coef_s1_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv')), data=matrix_coefs_csv_s1)
        call write_csv(filename=trim(StrStrStr(dir=trim(dir),str1='coef_s2_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv')), data=matrix_coefs_csv_s2)

        ! TEST IF SAVING LOSES DATA
        if (verbose) then !mute if not debugging
            nlines = size(matrix_csv_v1, 1)
            call read_csv(filename=StrStrStr(dir=trim(dir),str1='v1_data_'//mom2loop//'_', cnt=cnt_mom,str2='.csv'), data=matrix_csv_v1, nlines=nlines)

            v1_test = RESHAPE(source=(matrix_csv_v1), shape=shape(v1_test))
            meanval_v1 = maxval(abs(v1_test - v1))
            print *, '\033[34m MAX DISCREPANCY RELATIVE TO SAVED DATA \033[0m', meanval_v1

        end if

    end subroutine

end module
