module initialize_and_save_mod_EZ_k
    use my_kinds_mod, only: wp
    use globals
    implicit none
contains
! IF CRASHES HERE COULD BE BECAUSE OF MEMORY; DID YOU RUN THE LINECOMMAND ulimited -s unlimited ? 
    subroutine initialvalues_EZ_k(dir,cntry)
        use omp_lib
        use my_I_O
        use globals
        use read_csv_mod
        use my_kinds_mod
        use optimization_interface_mod
        implicit none

        integer :: ik1, ikp, iA1, iA2, is1, is2
        real(kind=wp):: kappa,lb,ub,parin0(7)
        real(kind=wp) :: Omega1, Omega2, C1, C2, mu, OmegaSplus,pf,ph,res,pred,L1,L2,Q,Y1
        character(len=5), intent(in) :: cntry
        character(len=100) :: def_filename, filename,dir
        logical :: file_exists, def_file_exists
        pf_sols = 1.0e0_wp
        if (doinitval) then
            !$omp PARALLEL DO  private(Omega1,Omega2,OmegaSplus,C1,C2,L1,L2,Y1,Q,res,pred,ph,pf,kappa,lb,ub,parin0) default(shared) collapse(6)

            do ik1 = 0, NK1
                do ikp = 0, NK4
                    do iA1 = 0, NK2
                        do iA2 = 0, NK3
                            do is1 = 1, NS
                                do is2 = 1, NS
                                    kappa=kappaS(ikp)
                 Omega1 = beta*(A1(iA1))**(1e0_wp/theta - 1e0_wp)*(-A1(iA1)**rho*exp(eta1(is1))**(1e0_wp - 1e0_wp/gamma)/(1e0_wp - 1e0_wp/gamma))**(theta - 1)

                 Omega2 = beta*(A2(iA2))**(1e0_wp/theta - 1e0_wp)*(-A2(iA2)**rho*exp(eta2(is2))**(1e0_wp - 1e0_wp/gamma)/(1e0_wp - 1e0_wp/gamma))**(theta - 1)
                                    ! update Omega
                                    OmegaSplus = Omega2/Omega1*OmegaS(ik1)
                                    pf=1.0_wp;
                                    call core_model(pf=pf, kappa=kappa, OmegaS_in=OmegaSplus, C=C1, Cs=C2, L=L1, Ls=L2, &
                                    ph=ph, Q=Q, Y1=Y1, res=res, pred=pred, D=A1(iA1)**rho_A1*exp(eta1(is1)), Ds=A2(iA2)**rho_A2*exp(eta2(is2)))

                                    v1(ik1, iA1, iA2,ikp, is1, is2) = C1**(1e0_wp - 1e0_wp/gamma)/(1e0_wp - 1e0_wp/gamma)-chhi*L1**(1.0e0_wp+varphi)/(1.0e0_wp+varphi) 
                                    v2(ik1, iA1, iA2,ikp, is1, is2) =  C2**(1e0_wp - 1e0_wp/gamma)/(1e0_wp - 1e0_wp/gamma)-chhi*L2**(1.0e0_wp+varphi)/(1.0e0_wp+varphi)
                                    s1(ik1, iA1, iA2,ikp, is1, is2) = (C1**(-1e0_wp/gamma)*(C1 - ph*Y1))/(1.0_wp-beta)
                                    v1_aut(iA1, iA2, is1, is2)=v1(ik1, iA1, iA2,ikp, is1, is2)
                                    v2_aut(iA1, iA2, is1, is2)=v2(ik1, iA1, iA2,ikp, is1, is2)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            !$omp END PARALLEL DO        

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


    end subroutine initialvalues_EZ_k

! ##############################################################################################################

        subroutine save_results_EZ_k(dir,cnt_mom)
        use globals
        use my_I_O
        use read_csv_mod
        implicit none

        character(len=5), intent(in) :: cnt_mom
        character(len=100), intent(in) :: dir
        character(len=100) :: filename, dir2, direnq
        integer :: ik1, unitnumber1, unitnumber2, unitnumber3, nlines
        logical :: dir_exists
! print*,dir
! print*,cnt_mom
                dir2 = trim(dir)
                inquire (directory=trim(dir2), exist=dir_exists)
                if (dir_exists .neqv. .true.) then
                    direnq = 'mkdir '//trim(dir2)
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
                call write_csv(filename=filename, data=matrix_csv_v1)
                filename = trim(StrStrStr(dir=trim(dir), str1='v2_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv'))
                call write_csv(filename=filename, data=matrix_csv_v2)
                filename = trim(StrStrStr(dir=trim(dir), str1='s1_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv'))
                call write_csv(filename=filename, data=matrix_csv_s1)
                filename = trim(StrStrStr(dir=trim(dir), str1='pB_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv'))
                call write_csv(filename=filename, data=matrix_csv_pf)

                ! SAVE AUTARKY
                matrix_csv_aut_v1 = RESHAPE(source=v1_aut, &
                                            shape=shape(matrix_csv_aut_v1))
                matrix_csv_aut_v2 = RESHAPE(source=v2_aut, &
                                            shape=shape(matrix_csv_aut_v2))

                matrix_csv_aut_pf = RESHAPE(source=pf_sols_aut, &
                                            shape=shape(matrix_csv_aut_pf))
                filename = trim(StrStrStr(dir=trim(dir), str1='v1_aut_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv'))
                call write_csv(filename=filename, data=matrix_csv_aut_v1)
                filename = trim(StrStrStr(dir=trim(dir), str1='v2_aut_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv'))
                call write_csv(filename=filename, data=matrix_csv_aut_v2)
                filename = trim(StrStrStr(dir=trim(dir), str1='pB_aut_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv'))
                call write_csv(filename=filename, data=matrix_csv_aut_pf)
                ! save coefficients
                matrix_coefs_csv_v1 = RESHAPE(source=coeff_v1, &
                                              shape=shape(matrix_coefs_csv_v1))
                matrix_coefs_csv_v2 = RESHAPE(source=coeff_v2, &
                                              shape=shape(matrix_coefs_csv_v2))

                matrix_coefs_csv_s1 = RESHAPE(source=coeff_s1, &
                                              shape=shape(matrix_coefs_csv_s1))
                filename = trim(StrStrStr(dir=trim(dir), str1='coef_v1_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv'))
                call write_csv(filename=filename, data=matrix_coefs_csv_v1)
                filename = trim(StrStrStr(dir=trim(dir), str1='coef_v2_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv'))
                call write_csv(filename=filename, data=matrix_coefs_csv_v2)
                filename = trim(StrStrStr(dir=trim(dir), str1='coef_s1_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv'))
                call write_csv(filename=filename, data=matrix_coefs_csv_s1)

                ! TEST IF SAVING LOSES DATA
                if (verbose) then !mute if not debugging
                    nlines = size(matrix_csv_v1, 1)
                    filename = StrStrStr(dir=trim(dir), str1='v1_data_'//mom2loop//'_', cnt=cnt_mom, str2='.csv')
                    call read_csv(filename=filename, data=matrix_csv_v1, nlines=nlines)

                    v1_test = RESHAPE(source=(matrix_csv_v1), shape=shape(v1_test))
                    meanval_v1 = maxval(abs(v1_test - v1))
                    print *, '\033[34m MAX DISCREPANCY RELATIVE TO SAVED DATA \033[0m', meanval_v1

                end if

            end subroutine

end module
