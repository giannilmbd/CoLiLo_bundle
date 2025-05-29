module initialize_and_save_mod
    use my_kinds_mod, only: wp
    implicit none
    contains

    subroutine initialvalues(kappa)
        use omp_lib
        use my_I_O
        use globals
        use read_csv_mod
        implicit none
    
        integer :: ik1,iA1,iA2,is1,is2
        real(kind=wp), intent(in) :: kappa
        real(kind=wp) :: Omega1,Omega2,C1,C2,mu,Yworld,OmegaSplus
        if (doinitval) then
            !$omp PARALLEL DO  private(Omega1,Omega2,OmegaSplus,YWorld,mu,C1,C2) default(shared) collapse(5)
    
            do ik1 = 0, NK1
                do iA1 = 0, NK2
                    do iA2 = 0, NK3
                        do is1 = 1, NS
                            do is2 = 1, NS
                 Omega1 = beta*(A1(iA1))**(1e0_wp/theta - 1e0_wp)*(-A1(iA1)**rho*exp(eta1(is1))**(1e0_wp - 1e0_wp/gamma)/(1e0_wp - 1e0_wp/gamma))**(theta - 1)
    
                 Omega2 = beta*(A2(iA2))**(1e0_wp/theta - 1e0_wp)*(-A2(iA2)**rho*exp(eta2(is2))**(1e0_wp - 1e0_wp/gamma)/(1e0_wp - 1e0_wp/gamma))**(theta - 1)
                                ! update Omega
                                OmegaSplus = Omega2/Omega1*OmegaS(ik1)
                                YWorld = n*A1(iA1)**rho*exp(eta1(is1)) + (1e0_wp - n)*A2(iA2)**rho*exp(eta2(is2))
                                mu = 1e0_wp/((kappa*OmegaSplus)**(gamma)*(1e0_wp - n) + n)
                                C1 = mu*YWorld
                                C2 = (1e0_wp/(1e0_wp - n) - mu*n/(1e0_wp - n))*YWorld
                                ! tmp_plot(ik1, iA1, iA2, is1, is2, 1) = C1**(1e0_wp - 1e0_wp/gamma)/(1e0_wp - 1e0_wp/gamma)
                                ! tmp_plot(ik1, iA1, iA2, is1, is2, 2) = OmegaS(ik1)
    
                                ! print*,'C1plot',C1plot
                                !G1
    
                                ! future X state
                                v1(ik1, iA1, iA2, is1, is2) = C1**(1e0_wp - 1e0_wp/gamma)/(1e0_wp - 1e0_wp/gamma)  ! check sidn Rud&Swans
                                v2(ik1, iA1, iA2, is1, is2) = C2**(1e0_wp - 1e0_wp/gamma)/(1e0_wp - 1e0_wp/gamma)
                                s1wide(ik1, iA1, iA2, is1, is2) = Omega1*(C1**(-1e0_wp/gamma)*(C1 - A1(iA1)))
                            end do
                        end do
                    end do
                end do
            end do
    !$omp END PARALLEL DO
    
        else
            call read_csv(filename=trim(StrNumStr('v1_data_'//mom2loop,1,'.csv')),data=matrix_csv_1)
            call read_csv(filename=trim(StrNumStr('v2_data_'//mom2loop,1,'.csv')),data=matrix_csv_2)
            call read_csv(filename=trim(StrNumStr('s1_data_'//mom2loop,1,'.csv')),data=matrix_csv_3)
            ! OPEN (unit=10, access="sequential", action="read", &
            ! file=trim(StrNumStr('v1_data_'//mom2loop,1,'.csv')) , form="formatted", status="old")
            ! OPEN (unit=11, access="sequential", action="read", &
            ! file=trim(StrNumStr('v2_data_'//mom2loop,1,'.csv')), form="formatted", status="old")
            ! OPEN (unit=12, access="sequential", action="read", &
            ! file=trim(StrNumStr('s1_data_'//mom2loop,1,'.csv')), form="formatted", status="old")
    
            ! do ik1 = 1, (NK2 + 1)*(NK3 + 1)*(NK1 + 1)*NS
            !     READ (10, *) matrix_csv_1(ik1, :)
            !     READ (11, *) matrix_csv_2(ik1, :)
            !     READ (12, *) matrix_csv_3(ik1, :)
    
            ! end do
            ! close (10)
            ! close (11)
            ! close (12)
            v1 = RESHAPE(source=(matrix_csv_1), shape=shape(v1))
            v2 = RESHAPE(source=(matrix_csv_2), shape=shape(v2))
            s1wide = RESHAPE(source=(matrix_csv_3), shape=shape(s1wide))
        end if
    ! print*,size(v1,1)
        ! tmp_plot2(:, 2) = RESHAPE(tmp_plot(:, :, :, :, :, 2), (/((1 + NK1)*(1 + NK2)*(1 + NK3)*NS**2)/))
        ! tmp_plot2(:, 1) = RESHAPE(tmp_plot(:, :, :, :, :, 1), (/((1 + NK1)*(1 + NK2)*(1 + NK3)*NS**2)/))
        ! call plot(tmp_plot2(:, 2), tmp_plot2(:, 1), noline=.true.)
        ! call execplot(xlabel='State', ylabel='C1', title='CONSUMPTION DISTRIBUTION (steady state)')
        v1_aut = v1(1, :, :, :, :); 
        v2_aut = v2(1, :, :, :, :);
        
        
    end subroutine initialvalues

    subroutine save_results(cnt_mom)
        use globals
        use my_I_O
        use read_csv_mod
        implicit none
    
        integer, intent(in) :: cnt_mom
        integer :: ik1,unitnumber1,unitnumber2,unitnumber3,nlines
        ! 101 format(1x, *(g0, ", "))
    ! SAVE THE SOLUTION TO FILE values and coefficients
        matrix_csv_1 = RESHAPE(source=v1, &
                               shape=shape(matrix_csv_1))
        matrix_csv_2 = RESHAPE(source=v2, &
                               shape=shape(matrix_csv_2))
        matrix_csv_3 = RESHAPE(source=s1wide, &
                               shape=shape(matrix_csv_3))
    
         call write_csv(filename=trim(StrNumStr('v1_data_'//mom2loop, cnt_mom, '.csv')),data=matrix_csv_1)
         call write_csv(filename=trim(StrNumStr('v2_data_'//mom2loop, cnt_mom, '.csv')),data=matrix_csv_2)
         call write_csv(filename=trim(StrNumStr('s1_data_'//mom2loop, cnt_mom, '.csv')),data=matrix_csv_3)

         ! SAVE AUTARKY
         matrix_csv_aut_1 = RESHAPE(source=v1_aut, &
         shape=shape(matrix_csv_aut_1))
         matrix_csv_aut_1 = RESHAPE(source=v2_aut, &
         shape=shape(matrix_csv_aut_2))


call write_csv(filename=trim(StrNumStr('v1_aut_data_'//mom2loop, cnt_mom, '.csv')),data=matrix_csv_aut_1)
call write_csv(filename=trim(StrNumStr('v2_aut_data_'//mom2loop, cnt_mom, '.csv')),data=matrix_csv_aut_2)
         
        ! save coefficients
    matrix_coefs_csv_1 = RESHAPE(source=coeff_v1, &
        shape=shape(matrix_coefs_csv_1))
    matrix_coefs_csv_2 = RESHAPE(source=coeff_v2, &
        shape=shape(matrix_coefs_csv_2))
    matrix_coefs_csv_3 = RESHAPE(source=coeff_s1, &
        shape=shape(matrix_coefs_csv_3))
        call write_csv(filename=trim(StrNumStr('coef_v1_data_'//mom2loop, cnt_mom, '.csv')),data=matrix_coefs_csv_1)
        call write_csv(filename=trim(StrNumStr('coef_v2_data_'//mom2loop, cnt_mom, '.csv')),data=matrix_coefs_csv_2)
        call write_csv(filename=trim(StrNumStr('coef_s1_data_'//mom2loop, cnt_mom, '.csv')),data=matrix_coefs_csv_3)
    ! OPEN (unit=10, access="sequential", action="write", &
    ! status="replace", &
    ! file=trim(StrNumStr('coef_v1_data_'//mom2loop, cnt_mom, '.csv')) &
    ! , form="formatted")
    ! OPEN (unit=11, access="sequential", action="write", &
    ! status="replace", &
    ! file=trim(StrNumStr('coef_v2_data_'//mom2loop, cnt_mom, '.csv')) &
    ! , form="formatted")
    ! OPEN (unit=12, access="sequential", action="write", &
    ! status="replace", &
    ! file=trim(StrNumStr('coef_s1_data_'//mom2loop, cnt_mom, '.csv')) &
    ! , form="formatted")
    ! ! Loop across rows
    ! do ik1 = 1,size(matrix_coefs_csv_1,1)
    ! WRITE (10, 101) matrix_coefs_csv_1(ik1, :)
    ! WRITE (11, 101) matrix_coefs_csv_2(ik1, :)
    ! WRITE (12, 101) matrix_coefs_csv_3(ik1, :)
    ! end do
    ! ! Close connection
    ! CLOSE (10)
    ! CLOSE (11)
    ! CLOSE (12)

    ! TEST IF SAVING LOSES DATA
    if(.true.)then !mute if not debugging
        nlines=size(matrix_csv_1,1)
        call read_csv(filename=StrNumStr('v1_data_'//mom2loop, cnt_mom, '.csv'), data=matrix_csv_1, nlines=nlines)
    
        v1_test = RESHAPE(source=(matrix_csv_1), shape=shape(v1_test))
        meanval_v1=maxval(abs(v1_test-v1))
        print*,'\033[34m MAX DISCREPANCY RELATIVE TO SAVED DATA \033[0m', meanval_v1

    endif 
    
    end subroutine


end module