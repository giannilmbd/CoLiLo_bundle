program main
   use globals
   use my_toolbox, only: toc, tic, init_random_seed, grid_Cons_Equi, fzero, simulate_AR, settol_root
   use omp_lib
   use discretize
   use mymoments
   use csv_module
   use read_csv_mod
   use linespaced_mod
   use my_I_O
   use initialize_and_save_mod_CRRA_2c
   use optimization_interface_mod
   use minpack_module, only: hybrd1, dpmpar, enorm
   use find_kappa_mod2, only: find_zero_arrow_bisection
   use pcu_CRRA_module, only: pcu_CRRA
   use MarkovMoments_mod, only: compute_moment
   USE, INTRINSIC :: IEEE_ARITHMETIC

   implicit none

   real(kind=wp) :: L1, L2, p1, Q, Y1, res, pred
   real(kind=wp):: xin1, xin2, C1, C2, Omega1, Omega2, OmegaSplus, welf_mom(n_moms_, 4)
   real(kind=wp), dimension(3) :: x_both, resids, XGUESS, FVEC
   real(kind=wp):: kappa, bis_a, bis_b, bis_x, resid_a, resid_b, resid, pops(actual_cntry, 2),moment1,moment2,moment3,moment4
   integer :: ires, nevals, indx_smpl, rc
   logical :: itexists
   integer :: cnt, fu, szrng, solver, slv, num_threads, cnt_mom, unitnumber1, unitnumber2, unitnumber3, unitnumber4
   integer, parameter :: nroots = 1
   real(kind=8):: EPS, ERRABS, ERRREL, kappa_roots(nroots)
   integer(kind=wp):: opt
   character(len=256), dimension(20) :: diagnostics
   real(kind=wp):: lb(3), minf, mu, YWorld, FNORM, parin0(100), moms_val_it
   integer :: IPARAM(6), countsim(TT)
   real(kind=wp) :: RPARAM(5), PDV, sigma_eps_tot, lbk(1), ubk(1), dummykappa(1),rngshk
   character*32:: date
   character*5:: cntry, def_cntry, cntry_b
   character*10 :: time
   character*100 :: dir
   character*30 :: starting
   integer, parameter :: size_welf = 18
   character*100, dimension(size_welf) :: header
   character*100, dimension(1, n_cols_ar) :: header_ar
   real(kind=wp), dimension(n_moms_, size_welf) :: welfout
   
   character(len=30),dimension(:),allocatable:: header_share
   character(len=30),dimension(:),allocatable :: countries_share2
   character(len=5),dimension(N_countries) :: countries_share
   real(kind=wp),dimension(:),allocatable :: cons_imp_share
   real(kind=wp):: pr1, pr2, pr3, var1, var2, var3, mu1, mu2, mu3, m4
   integer ::  ik1, ik2, ik3, is1, is2
   real(kind=wp):: kappa_max_iter, exit_crit, resid_change,pcu_c1_au
   character::tmpfmt
   character*100:: file_data, tmpmoms
   external F_kappa
   type(csv_file) :: f
   character(len=30), dimension(:), allocatable :: header_new
   real(wp), dimension(:), allocatable :: ar_stats2
   character(len=100), dimension(:), allocatable :: name_countries2
   logical :: status_ok
   integer :: cntry_1
   integer, dimension(:), allocatable :: itypes
   real(kind=wp) :: v1_aut_init, v2_aut_init

   write (*, *) ANSI_BLUE//'========================================================'
   write (*, *) ANSI_BLUE//'Allocating matrices'
   write (*, *) ANSI_BLUE//'========================================================'//ANSI_RESET
   call allocate_globals()
   write (*, *) ANSI_YELLOW//'========================================================'
   write (*, *) ANSI_YELLOW//'Allocating matrices done'
   write (*, *) ANSI_YELLOW//'========================================================'//ANSI_RESET

   call date_and_time(date, time)
   starting = 'On '//date(:4)//'.'//date(5:6)//'.'//date(7:8)//' at '//time(1:2)//':'//time(3:4)//':'//time(5:)

   IF (COMMAND_ARGUMENT_COUNT() .NE. 2) THEN
      WRITE (*, '(a)') ANSI_RED//'ERROR, COMMAND-LINE ARGUMENTS REQUIRED: (v,s,k) and number of moments, STOPPING'//ANSI_RESET
      STOP
   END IF
   call get_command_argument(number=1, value=mom2loop, status=rc)
   call get_command_argument(number=2, value=tmpmoms, status=rc)
   read (tmpmoms, *) n_moms
   n_moms = actual_cntry
   WRITE (*, '(a,i5)') ANSI_YELLOW//'NUMBER OF MOMENTS TO CONSIDER'//ANSI_RESET, n_moms

   call init_random_seed(fixed)
   call settol_root(1e-7_wp)

   OmegaS = linespaced_values(NK4 + 1, LB_all(3), UB_all(3))
   OmegaS(INT((NK4+1)/2.0))=1.0_wp
   pos_zeroOm = (INT((NK4+1)/2.0))

   write (*, '(A87)') adjustl('[91m ############################################################################ [0m')
   write (*, '(A87)') adjustr('[91m ####               Running session with parameters                [0m')
   write (*, '(A55,A30,A5)') adjustl('[91m ####               Started                                         '), starting, ' [0m'
   write (*, '(A42,I4,A5)') adjustl('[91m ####               Nbr moments:                                    '), n_moms, ' [0m'
   write (*, '(A42,A1,A5)') adjustl('[91m ####               Type moment:                                    '), mom2loop, ' [0m'
   write (*, '(A55,L1,A5)') adjustl('[91m ####               Ad-hoc initial conditions:                      '), doinitval, ' [0m'
   write (*, '(A55,A20,A5)') adjustl('[91m ####               Treatment of Euler Equations:                   '), inner_loop, ' [0m'
   write (*, '(A87)') adjustl('[91m ####               Region:                   '//region//' [0m')
   write (*, '(A87)') adjustl('[91m ############################################################################ [0m')

   call tic()

   call f%read('../../../Rcode/ar_data_pwt.csv', header_row=1, status_ok=status_ok)

   call f%get_header(header_new, status_ok)
   call f%variable_types(itypes, status_ok)

   call f%get(1, name_countries2, status_ok)
   call f%get(2, ar_stats2, status_ok)
   
   print *, "Size of ar_stats2:", size(ar_stats2)
   print *, "Expected size:", N_countries
   
   if (size(ar_stats2) /= N_countries) then
      print *, "ERROR: Size mismatch in ar_stats2"
      print *, "Expected:", N_countries, "Got:", size(ar_stats2)
      stop
   end if
   
   ar_stats(:, 1) = ar_stats2
   do cnt = 2, n_cols_ar
      call f%get(cnt, ar_stats2, status_ok)
      if (size(ar_stats2) /= N_countries) then
         print *, "ERROR: Size mismatch in ar_stats2 at column", cnt
         print *, "Expected:", N_countries, "Got:", size(ar_stats2)
         stop
      end if
      ar_stats(:, cnt - 1) = ar_stats2
   end do

   do cnt = 1, N_countries
      names_countries(cnt) = trim(name_countries2(cnt))
   end do

   call f%destroy()

   do cnt = 1, actual_cntry
      pos_subsample(cnt) = FINDLOC(names_countries, subsample(cnt), 1)
   end do

   call f%read('../../../Rcode/ConsumptionImportShares.csv',header_row=1,status_ok=status_ok)

   call f%get_header(header_share,status_ok)
   call f%variable_types(itypes,status_ok)

   call f%get(1,countries_share2,status_ok)
   call f%get(2,cons_imp_share,status_ok)
   do cnt = 1, N_countries
      countries_share(cnt)=trim(countries_share2(cnt))
   end do
   call f%destroy()

   tmpmoms = '('//trim(tmpmoms)//'e12.5)'
   def_cntry =trim(names_countries(12))

   inquire (file=StrNumStr(str1='../../../Calibrate_moments_py/Pmatrices_iid_', cnt=NS, str2='.csv'), exist=status_ok)
   if (status_ok .neqv. .true.) then
      print*, ANSI_RED//'NO FILE ==> '//StrNumStr(str1='../../Calibrate_moments_py/Pmatrices_iid_', cnt=NS, str2='.csv')//' STOPPING'//ANSI_RESET
      stop
   end if
   call read_csv(filename=StrNumStr(str1='../../../Calibrate_moments_py/Pmatrices_iid_', cnt=NS, str2='.csv'), data=allP)

   call read_csv(filename=StrNumStr(str1='../../../Calibrate_moments_py/Xmatrices_iid_', cnt=NS, str2='.csv'), data=allX)

   dir = adjustl(trim('../csvfiles'//trim(region)))

   loop_moments: do cnt_mom = 1, n_moms
      write (*, '(A40,1x,i5,A40,1x,A40)') ANSI_BLUE//' moment number  '//ANSI_RESET, cnt_mom, ANSI_YELLOW//'country'//ANSI_RESET ,ANSI_BACKGROUND//ANSI_YELLOW//trim(names_countries(pos_subsample(cnt_mom)))//ANSI_RESET

      if (region == '_D') then
         base_cntry = FINDLOC(names_countries, 'DEU', 1)
      elseif (region == '_U') then
         base_cntry = FINDLOC(names_countries, 'USA', 1)
      else
         base_cntry = FINDLOC(names_countries, trim(names_countries(pos_subsample(cnt_mom)))//region, 1)
      end if
      write (*, *) ANSI_RED//ANSI_BACKGROUND//'BASE COUNTRY '//ANSI_RESET, names_countries(base_cntry)
      pi2 = transpose(reshape(source=allP(base_cntry, :), shape=shape(pi2)))
      eta2 = reshape(source=allX(base_cntry, :), shape=shape(eta2))
      pop2 = ar_stats(base_cntry, n_cols_ar - 1)
      rho_A2 = ar_stats(base_cntry,1)

      pos_zero2 = minloc(abs(eta2), 1)

      pi1 = transpose(reshape(source=allP(pos_subsample(cnt_mom), :), shape=shape(pi1)))
      eta1 = reshape(source=allX(pos_subsample(cnt_mom), :), shape=shape(eta1))

      pos_zero1 = minloc(abs(eta1), 1)

      mom_val_1(cnt_mom, 1) =0

      call compute_moment(pi1, eta1,2.0_wp,moment2)
      mom_val_1(cnt_mom, 2) = sqrt(moment2)
      call compute_moment(pi1, eta1,3.0_wp,moment3)
      mom_val_1(cnt_mom, 3) = moment3/moment2**(3.0/2.0)
      call compute_moment(pi1, eta1,4.0_wp,moment4)
      mom_val_1(cnt_mom, 4) = moment4/moment2**2.0

      mom_val_2(cnt_mom, 1) = 0
      call compute_moment(pi2, eta2,2.0_wp,moment2)
      mom_val_2(cnt_mom, 2) = sqrt(moment2)

      call compute_moment(pi2, eta2,3.0_wp,moment3)
      mom_val_2(cnt_mom, 3) = moment3/moment2**(3.0/2.0)

      call compute_moment(pi2, eta2,4.0_wp,moment4)
      mom_val_2(cnt_mom, 4) = moment4/moment2**2.0

      rho_A1 = ar_stats(pos_subsample(cnt_mom),1)
      rngshk=5.0_wp*mom_val_1(cnt_mom, 2)/sqrt(1.0_wp-rho_A1**2)

      LB_all(1)=exp(-rngshk)
      
      UB_all(1)=exp(rngshk)
      rngshk=5.0_wp*mom_val_2(cnt_mom, 2)/sqrt(1.0_wp-rho_A2**2)
      LB_all(2)=exp(-rngshk)
      UB_all(2)=exp(rngshk)
      print*,ANSI_YELLOW,"LB: ",LB_all,ANSI_RESET
      A1 = linespaced_values(NK2 + 1, LB_all(1), UB_all(1))
      A1(INT((NK2+1)/2.0))=1.0_wp
      A2 = linespaced_values(NK3 + 1, LB_all(2), UB_all(2))
      A2(INT((NK3+1)/2.0))=1.0_wp

      pos_zeroA1 = INT((NK2+1)/2.0)
      pos_zeroA2 = INT((NK3+1)/2.0)

      write (*, '(A5,A20,1x,4A20)') 'country|', 'type|', 'Mean|', 'Var|', 'Skewness|', 'kurtosis|'

      write (*, '(A5,A20,1x,4f20.10)') names_countries(pos_subsample(cnt_mom)),' Simulated', 0.0,mom_val_1(cnt_mom, 2)**2.0, mom_val_1(cnt_mom, 3), mom_val_1(cnt_mom, 4)

      write (*, '(A5,A20,1x,4f20.10)') names_countries(pos_subsample(cnt_mom)), 'Empirical', ar_stats(cnt_mom, 6), &
         ar_stats(cnt_mom, 7), &
         ar_stats(cnt_mom, 8), &
         ar_stats(cnt_mom, 9)
      write (*, '(A5,A20,1x,4f20.10)') names_countries(base_cntry),'Simulated',0.0,mom_val_2(cnt_mom, 2)**2.0, mom_val_2(cnt_mom, 3), mom_val_2(cnt_mom, 4)

      cntry = names_countries(pos_subsample(cnt_mom))
      pop1 = ar_stats(pos_subsample(cnt_mom), n_cols_ar - 1)
      n = pop1/(pop1 + pop2)
      nu = (1.0_wp-cons_imp_share(pos_subsample(cnt_mom))/100_wp)
      nu_b =(1.0_wp-cons_imp_share(base_cntry)/100_wp)

      pops(cnt_mom, 1) = pop1
      pops(cnt_mom, 2) = pop2
      price_down =1.0e-2_wp
      write (*, *) ANSI_BLUE//'========================================================'
      write (*, *) ANSI_BLUE//' Population size (A,B)', n, (1.0-n)
      write (*, *) ANSI_BLUE//' Domestic Goods Weight (A,B)', nu, nu_b
      write (*, *) ANSI_BLUE//' Lower Bound Price', price_down
      write (*, *) ANSI_BLUE//'========================================================'//ANSI_RESET

      write (*, *) ANSI_BLUE//'========================================================'
      write (*, *) ANSI_BLUE//'INITIALIZING'
      write (*, *) ANSI_BLUE//'========================================================'//ANSI_RESET
      call initialvalues_CRRA(dir, cntry,def_cntry)
      write (*, *) ANSI_BLUE//'========================================================'
      write (*, '(A,A)') ANSI_BLUE//'INITIALIZING DONE','✔'  
      write (*, *) ANSI_BLUE//'========================================================'//ANSI_RESET
      write (*,'(A,f15.6,A,f15.6,A)') ANSI_BLUE//'Welfare aut at ss at initialization: V1 = ',v1_aut(pos_zeroA1,pos_zeroA2,pos_zero1,pos_zero2),&
      "; V2 = ",v2_aut(pos_zeroA1,pos_zeroA2,pos_zero1,pos_zero2),ANSI_RESET
      ! print*, "Press enter to continue"
      ! read(*,*) cnt 
      write (*, *) ANSI_BLUE//'========================================================'
      write (*, *) ANSI_BLUE//'SOLVING AUTARKY'
      v1_aut_init=v1_aut(pos_zeroA1,pos_zeroA2,pos_zero1,pos_zero2)
      v2_aut_init=v2_aut(pos_zeroA1,pos_zeroA2,pos_zero1,pos_zero2)
      call value_autarky()
      write(*,*) ANSI_YELLOW, "v1 AUT INIT", v1_aut_init
      write(*,*) ANSI_YELLOW, "v2 AUT INIT", v2_aut_init

      write(*,*) ANSI_YELLOW, "v1 AUT", v1_aut(pos_zeroA1,pos_zeroA2,pos_zero1,pos_zero2)
      write(*,*) ANSI_YELLOW, "v2 AUT", v2_aut(pos_zeroA1,pos_zeroA2,pos_zero1,pos_zero2)

   end do loop_moments

   call deallocate_globals()
   write (*, *) ANSI_BLUE//'========================================================'
   write (*, '(A)') 'END OF SIMULATION ✔'
   
   call execute_command_line("date")
   write (*, '(A)') ANSI_BLUE//'========================================================'//ANSI_RESET

end program main

!gfortran -O3 -march=native -mtune=native -fopenmp  globals.f90 my_kinds_mod.f90 test_autarky.f90 -o test_autarky