
program save_coeffs_one_time_EZ_k

    use globals, disabled => foc2
    use toolbox
    use my_I_O
    use discretize
    use mymoments
    use read_csv_mod
    use linespaced_mod
    USE, INTRINSIC :: IEEE_ARITHMETIC
    
    implicit none

    real*8 :: xin1, xin2, C1, C2, Omega1, Omega2, OmegaSplus
    real*8, dimension(2) :: x_both, resids, XGUESS, FVEC
    real*8 :: kappa, bis_a, bis_b, bis_x, resid_a, resid_b
    integer ::  from, to, arglength, nunit, nlines, ires, nevals, indx_smpl, cnt2, rc
    character(len=100) :: dir, cntstr, tmpmoms
    character(len=:), allocatable :: degasy
    logical :: itexists
    integer :: cnt, fu, szrng, solver, slv, num_threads, cnt_mom


    integer ::  ik1, ik2, ik3, is1, is2

    character::tmpfmt
    character*100:: file_data
    


    IF (COMMAND_ARGUMENT_COUNT() .NE. 3) THEN
        WRITE (*, *) 'ERROR, COMMAND-LINE ARGUMENTS REQUIRED, STOPPING'
        STOP
    END IF
    call get_command_argument(number=1, value=mom2loop, status=rc)
    call get_command_argument(number=2, length=arglength, status=rc)
    call get_command_argument(number=3, value=dir, status=rc)
    allocate (character(arglength):: degasy)
    call get_command_argument(number=2, value=degasy, status=rc)
! call get_command_argument(number=3, value=kappainput, status=rc)

    read (degasy, *) cnt2

    cnt_mom = cnt2


    allocate(matrix_csv_1((NK2 + 1)*(NK3 + 1)*(NK1 + 1)*(NK4+1)*NS, NS))
    call read_csv(filename=trim(StrNumStr('v1_data_'//mom2loop, cnt_mom, '.csv')), data=matrix_csv_1)
    v1 = RESHAPE(source=(matrix_csv_1), shape=shape(v1))
    deallocate(matrix_csv_1)
    allocate(matrix_csv_2((NK2 + 1)*(NK3 + 1)*(NK1 + 1)*(NK4+1)*NS, NS))
    call read_csv(filename=trim(StrNumStr('v2_data_'//mom2loop, cnt_mom, '.csv')), data=matrix_csv_2)
    v1 = RESHAPE(source=(matrix_csv_2), shape=shape(v2))
    deallocate(matrix_csv_2)
    allocate(matrix_csv_3((NK2 + 1)*(NK3 + 1)*(NK1 + 1)*(NK4+1)*NS, NS))
    call read_csv(filename=trim(StrNumStr('s1_data_'//mom2loop, cnt_mom, '.csv')), data=matrix_csv_3)
    s1 = RESHAPE(source=(matrix_csv_3), shape=shape(s1))
    deallocate(matrix_csv_3)

!$omp PARALLEL  default(shared)
        !$omp DO collapse(2)
    do is1 = 1, NS
        do is2 = 1, NS

            call spline_interp(v1(:, :, :,:, is1, is2), coeff_v1(:, :, :, :,is1, is2))
            call spline_interp(v2(:, :, :,:, is1, is2), coeff_v2(:, :, :, :,is1, is2))
            call spline_interp(s1(:, :, :,:, is1, is2), coeff_s1(:, :, :,:, is1, is2))
            ! tmp_plot(:,:,:, is1, is2, 1) = v1(:,:,:,is1, is2)

            ! tmp_plot(:,:,:, is1, is2, 2) = spline_eval([OmegaS(ik1), exp(eta1(ik2)), exp(eta2(ik3))], &
            !       coeff_v1(:, :, :, is1, is2), LB_all, UB_all)
        end do
    end do
    !$omp END DO
    !$omp END PARALLEL
    
    call write_csv(filename=trim(StrNumStr('coef_v1_data_'//mom2loop, cnt_mom, '.csv')), data=RESHAPE(source=coeff_v1, shape=(/(NK2 + 3)*(NK3 + 3)*(NK1 + 3)*(NK4+3)*NS,NS/)))
    call write_csv(filename=trim(StrNumStr('coef_v2_data_'//mom2loop, cnt_mom, '.csv')), data=RESHAPE(source=coeff_v2, shape=(/(NK2 + 3)*(NK3 + 3)*(NK1 + 3)*(NK4+3)*NS,NS/)))

    call write_csv(filename=trim(StrNumStr('coef_s1_data_'//mom2loop, cnt_mom, '.csv')), data=RESHAPE(source=coeff_s1, shape=(/(NK2 + 3)*(NK3 + 3)*(NK1 + 3)*(NK4+3)*NS,NS/)))
   
    end program save_coeffs_one_time_EZ_k