! THIS VERSION READS THE EMPIRICAL - country-based - MOMENTS' PARAMETRIZATION FROM FILE AND LOOPS OVER THEM
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
   !  use info_window
   use MarkovMoments_mod, only: compute_moment
   ! USE ZREAL_INT
   ! USE WRRRN_INT
   USE, INTRINSIC :: IEEE_ARITHMETIC
   ! use ogpf
   ! use lapackMod

   implicit none

   ! real(kind=wp):: resids(2)

   real(kind=wp) :: L1, L2, p1, Q, Y1, res, pred
   real(kind=wp):: xin1, xin2, C1, C2, Omega1, Omega2, OmegaSplus, welf_mom(n_moms_, 4)!,parin(9)
   real(kind=wp), dimension(3) :: x_both, resids, XGUESS, FVEC
   real(kind=wp):: kappa, bis_a, bis_b, bis_x, resid_a, resid_b, resid, pops(actual_cntry, 2),moment1,moment2,moment3,moment4
   integer :: ires, nevals, indx_smpl, rc
   logical :: itexists
   integer :: cnt, fu, szrng, solver, slv, num_threads, cnt_mom, unitnumber1, unitnumber2, unitnumber3, unitnumber4
   integer, parameter :: nroots = 1
   real(kind=8):: EPS, ERRABS, ERRREL, kappa_roots(nroots)
   integer(kind=wp):: opt
    character(len=256), dimension(20) :: diagnostics
   ! TYPE(gpf):: gp
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
   ! real(kind=wp):: tmp_plot(0:NK1, 0:NK2, 0:NK3, NS, NS, 2), tmp_plot2((1 + NK1)*(1 + NK2)*(1 + NK3)*NS**2, 2)
    integer ::  ik1, ik2, ik3, is1, is2
    real(kind=wp):: kappa_max_iter, exit_crit, resid_change,pcu_c1_au
   character::tmpfmt
   character*100:: file_data, tmpmoms
    external F_kappa
   ! variables to read csv
   type(csv_file) :: f
   character(len=30), dimension(:), allocatable :: header_new
   real(wp), dimension(:), allocatable :: ar_stats2
   character(len=100), dimension(:), allocatable :: name_countries2
   logical :: status_ok
    integer :: cntry_1
   integer, dimension(:), allocatable :: itypes

   write (*, *) ANSI_BLUE//'========================================================'
   write (*, *) ANSI_BLUE//'Allocating matrices'
   write (*, *) ANSI_BLUE//'========================================================'//ANSI_RESET
call allocate_globals()
write (*, *) ANSI_YELLOW//'========================================================'
write (*, *) ANSI_YELLOW//'Allocating matrices done'
write (*, *) ANSI_YELLOW//'========================================================'//ANSI_RESET

!! ENQUIRE WHICH MOMENTS TO LOOP OVER
   call date_and_time(date, time)
   starting = 'On '//date(:4)//'.'//date(5:6)//'.'//date(7:8)//' at '//time(1:2)//':'//time(3:4)//':'//time(5:)
   ! write (*, '(a)') ANSI_RED//'CHOOSE MOMENTS TO LOOP OVER: v=VARIANCE; s=SKEWNESS; k=KURTOSIS '//ANSI_RESET
   ! !read (*, '(1A)') mom2loop
   ! in running the program must give number of countries and a character (eg v... legacy input to be removed in later versions)
   IF (COMMAND_ARGUMENT_COUNT() .NE. 2) THEN
      WRITE (*, '(a)') ANSI_RED//'ERROR, COMMAND-LINE ARGUMENTS REQUIRED: (v,s,k) and number of moments, STOPPING'//ANSI_RESET
      STOP
   END IF
   call get_command_argument(number=1, value=mom2loop, status=rc)
   call get_command_argument(number=2, value=tmpmoms, status=rc)
   read (tmpmoms, *) n_moms
   ! n_moms = min(n_moms, N_countries) ! override as use only hard-wired subsample
   n_moms = actual_cntry
   WRITE (*, '(a,i5)') ANSI_YELLOW//'NUMBER OF MOMENTS TO CONSIDER'//ANSI_RESET, n_moms
   !include 'nlopt.f'
100 format(*(E20.12, ", "))
101 format(1x, *(E20.12, ", "))
   call init_random_seed(fixed)
   call settol_root(1e-7_wp)
   ! INITIALIZE V1 and V2 (not clear why they give problems)
!    v1(:, :, :) = 1.0e0_wp
   !  v2(:, :, :) = 1.0e0_wp


   ! write(*,*) '[95m USED SAVED SOLUTION AS INITIAL CONDITIONS? (1/0)[0m '
   ! read(*,*) doinitial



   OmegaS = linespaced_values(NK4 + 1, LB_all(3), UB_all(3))
   OmegaS(INT((NK4+1)/2.0))=1.0_wp
   pos_zeroOm = (INT((NK4+1)/2.0))
   ! write(*,'(a,e20.12)') ANSI_BLUE//'KAPPA RANGE: '//ANSI_RESET, kappaS
   write (*, '(A87)') ANSI_RED//' ############################################################################ '//ANSI_RESET
   write (*, '(A87)') ANSI_RED//' ####               Running session with parameters                '//ANSI_RESET
   write (*, '(A55,A30,A5)') ANSI_RED//' ####               Started                                         ', starting, ' '//ANSI_RESET
   write (*, '(A42,I4,A5)') ANSI_RED//' ####               Nbr moments:                                    ', n_moms, ' '//ANSI_RESET
   write (*, '(A42,A1,A5)') ANSI_RED//' ####               Type moment:                                    ', mom2loop, ' '//ANSI_RESET
   write (*, '(A55,L1,A5)') ANSI_RED//' ####               Ad-hoc initial conditions:                      ', doinitval, ' '//ANSI_RESET
   write (*, '(A55,A20,A5)') ANSI_RED//' ####               Treatment of Euler Equations:                   ', inner_loop, ' '//ANSI_RESET
   write (*, '(A87)') ANSI_RED//' ####               Region:                   '//region//ANSI_RESET
   write (*, '(A87)') ANSI_RED//' ############################################################################ '//ANSI_RESET
!################################################ INITIALIZE POLICY FUNCTIONS ##########################################

   !print*,'IN MAIN INITIAL S1', s1
   call tic()



! LOAD MARKOV PROCESS
!# LOAD AR COEFFS
   call f%read('../../Rcode/ar_data_pwt.csv', header_row=1, status_ok=status_ok)

   ! get the header and type info
   call f%get_header(header_new, status_ok)
   call f%variable_types(itypes, status_ok)

   ! get some data
   call f%get(1, name_countries2, status_ok)
   call f%get(2, ar_stats2, status_ok)
   
   ! Debug print
   print *, "Size of ar_stats2:", size(ar_stats2)
   print *, "Expected size:", N_countries
   
   ! Check if sizes match
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
      ar_stats(:, cnt - 1) = ar_stats2 ! NOTE THAT HERE IS -1 SO THE SIZE WILL BE CONSISTENT
   end do

   do cnt = 1, N_countries
      names_countries(cnt) = trim(name_countries2(cnt))
   end do
   ! names_countries = trim(name_countries2)
   ! destroy the file
   call f%destroy()

! find position of subsample
   do cnt = 1, actual_cntry
      pos_subsample(cnt) = FINDLOC(names_countries, subsample(cnt), 1)

   end do

   ! LOAD TRADE OPENNESS
   call f%read('../../Rcode/ConsumptionImportShares.csv',header_row=1,status_ok=status_ok)

   ! get the header and type info
   call f%get_header(header_share,status_ok)
   call f%variable_types(itypes,status_ok)

   ! get some data
   call f%get(1,countries_share2,status_ok)
   call f%get(2,cons_imp_share,status_ok)
   do cnt = 1, N_countries
   countries_share(cnt)=trim(countries_share2(cnt))
   end do
   call f%destroy()

   tmpmoms = '('//trim(tmpmoms)//'e12.5)'
   def_cntry =trim(names_countries(12))
! LOAD MARKOV P and X
   inquire (file=StrNumStr(str1='../../Calibrate_moments_py/Pmatrices_iid_', cnt=NS, str2='.csv'), exist=status_ok)
   if (status_ok .neqv. .true.) then
      print*, ANSI_RED//'NO FILE ==> '//StrNumStr(str1='../../Calibrate_moments_py/Pmatrices_iid_', cnt=NS, str2='.csv')//' STOPPING'//ANSI_RESET
      stop
   end if
   call read_csv(filename=StrNumStr(str1='../../Calibrate_moments_py/Pmatrices_iid_', cnt=NS, str2='.csv'), data=allP)

   call read_csv(filename=StrNumStr(str1='../../Calibrate_moments_py/Xmatrices_iid_', cnt=NS, str2='.csv'), data=allX)

   ! call simulate(Cpimx2, zindex2) ! FROM BOOK
   ! ! totindex=[INT(zindex),NS+1-INT(zindex)]
   ! totindex2 = INT(zindex2)
   ! eta_simul2 = eta2(totindex2)
   ! sigma_eps_tot = p1_ar*(m1**2.0e0_wp + inputpar1(2)) + p2_ar*(m2**2.0e0_wp + inputpar1(4)) + p3_ar*(m3**2.0e0_wp + inputpar1(6))
   ! print *, 'sigma total', sigma_eps_tot
   ! use same base for distribution across countries

   !######################################### LOOP MOMENTS #####################################

   dir = adjustl(trim('../csvfiles'//trim(region)))


   loop_moments: do cnt_mom = 1, n_moms
      ! mom_range(p1,p2,mu2,v1,v2)
      write (*, '(A40,1x,i5,A40,1x,A40)') ANSI_BLUE//' moment number  '//ANSI_RESET, cnt_mom, ANSI_YELLOW//'country'//ANSI_RESET ,ANSI_BACKGROUND//ANSI_YELLOW//trim(names_countries(pos_subsample(cnt_mom)))//ANSI_RESET
      ! FOR GERMANY OR USA BETTER TO HARD WIRE THE DISTRIBUTION OF THESE COUNTRIES RATHER THAN RELYING ON NEW MATCHING ()
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
      pop2 = ar_stats(base_cntry, n_cols_ar - 1) !*2
      rho_A2 = ar_stats(base_cntry,1)
      ! print*,ANSI_RED//'NB: CHANGED POPULATION 2 To TWICE WHAT IT IS FOR 2c vs 3c comparison'//ANSI_RESET

      pos_zero2 = minloc(abs(eta2), 1)
      ! call simulate_AR(pi2, totindex2, .true.)
      ! eta_simul2 = eta2(totindex2)

      pi1 = transpose(reshape(source=allP(pos_subsample(cnt_mom), :), shape=shape(pi1)))
      eta1 = reshape(source=allX(pos_subsample(cnt_mom), :), shape=shape(eta1))

      pos_zero1 = minloc(abs(eta1), 1)
      ! call simulate_AR(pi1, totindex1, .true.)
      ! eta_simul1 = eta1(totindex1)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! HOME
      mom_val_1(cnt_mom, 1) =0

      call compute_moment(pi1, eta1,2.0_wp,moment2)
      mom_val_1(cnt_mom, 2) = sqrt(moment2)
      call compute_moment(pi1, eta1,3.0_wp,moment3)
      mom_val_1(cnt_mom, 3) = moment3/moment2**(3.0/2.0)
      call compute_moment(pi1, eta1,4.0_wp,moment4)
      mom_val_1(cnt_mom, 4) = moment4/moment2**2.0
! FOREIGN
      mom_val_2(cnt_mom, 1) = 0
      call compute_moment(pi2, eta2,2.0_wp,moment2)
      mom_val_2(cnt_mom, 2) = sqrt(moment2)

      call compute_moment(pi2, eta2,3.0_wp,moment3)
      mom_val_2(cnt_mom, 3) = moment3/moment2**(3.0/2.0)

      call compute_moment(pi2, eta2,4.0_wp,moment4)
      mom_val_2(cnt_mom, 4) = moment4/moment2**2.0


      rho_A1 = ar_stats(pos_subsample(cnt_mom),1)
      ! bounds should be commensurate to the shocks eg 5 stdev
      rngshk=5.0_wp*mom_val_1(cnt_mom, 2)/sqrt(1.0_wp-rho_A1**2)

      LB_all(1)=exp(-rngshk)
      
      UB_all(1)=exp(rngshk)
      rngshk=5.0_wp*mom_val_2(cnt_mom, 2)/sqrt(1.0_wp-rho_A2**2)
      LB_all(2)=exp(-rngshk)
      UB_all(2)=exp(rngshk)
print*,ANSI_YELLOW,"LB: ",LB_all,ANSI_RESET
      A1 = linespaced_values(NK2 + 1, LB_all(1), UB_all(1))
      A1(INT((NK2+1)/2.0))=1.0_wp
      ! print*, "A1",A1
      A2 = linespaced_values(NK3 + 1, LB_all(2), UB_all(2))
      A2(INT((NK3+1)/2.0))=1.0_wp


! MAKE IT SYMMETRIC
      
! A2=A1
! UB_all(3)=UB_all(2)
! LB_all(3)=LB_all(2)
! pi2=pi1
! eta2=eta1      

!  print*,A1
!  print*,A2
!  print*,OmegaS
! stop
!    call grid_Cons_Equi(OmegaS, Omega_l, Omega_u)
!   call grid_Cons_Equi(A1, A1_l, A1_u)
!   call grid_Cons_Equi(A2, A2_l, A2_u)

      pos_zeroA1 = INT((NK2+1)/2.0)! minloc(abs(A1 - 1.0e0_wp), 1) - 1
      pos_zeroA2 = INT((NK3+1)/2.0) !minloc(abs(A2 - 1.0e0_wp), 1) - 1


      write (*, '(A5,A20,1x,4A20)') 'country|', 'type|', 'Mean|', 'Var|', 'Skewness|', 'kurtosis|'

      write (*, '(A5,A20,1x,4f20.10)') names_countries(pos_subsample(cnt_mom)),' Simulated', 0.0,mom_val_1(cnt_mom, 2)**2.0, mom_val_1(cnt_mom, 3), mom_val_1(cnt_mom, 4)

      write (*, '(A5,A20,1x,4f20.10)') names_countries(pos_subsample(cnt_mom)), 'Empirical', ar_stats(cnt_mom, 6), &
         ar_stats(cnt_mom, 7), &
         ar_stats(cnt_mom, 8), &
         ar_stats(cnt_mom, 9)
      write (*, '(A5,A20,1x,4f20.10)') names_countries(base_cntry),'Simulated',0.0,mom_val_2(cnt_mom, 2)**2.0, mom_val_2(cnt_mom, 3), mom_val_2(cnt_mom, 4)



      cntry = names_countries(pos_subsample(cnt_mom))
      ! Determine population size
      pop1 = ar_stats(pos_subsample(cnt_mom), n_cols_ar - 1)
      n = pop1/(pop1 + pop2)
      nu = (1.0_wp-cons_imp_share(pos_subsample(cnt_mom))/100_wp)! (1.0_wp - (1.0_wp - n)*lamb)
      nu_b =(1.0_wp-cons_imp_share(base_cntry)/100_wp)!(1.0_wp - ( n)*lamb)!

      pops(cnt_mom, 1) = pop1
      pops(cnt_mom, 2) = pop2
      price_down =1.0e-2_wp! (1.0e0_wp - nu)**(1/(trade_elast - 1.0e0_wp)) + (1.0e-8_wp) ! ditto
      write (*, *) ANSI_BLUE//'========================================================'
      write (*, *) ANSI_BLUE//' Population size (A,B)', n, (1.0-n)
      write (*, *) ANSI_BLUE//' Domestic Goods Weight (A,B)', nu, nu_b
      write (*, *) ANSI_BLUE//' Lower Bound Price', price_down
      write (*, *) ANSI_BLUE//'========================================================'//ANSI_RESET

      ! the default counry  is set in globals: used only if the initial conditions for cntry are not found. The results will be for cntry nevertheless
      write (*, *) ANSI_BLUE//'========================================================'
      write (*, *) ANSI_BLUE//'INITIALIZING'
      write (*, *) ANSI_BLUE//'========================================================'//ANSI_RESET
      call initialvalues_CRRA(dir, cntry,def_cntry)
      write (*, *) ANSI_BLUE//'========================================================'
      write (*, '(A,A)') ANSI_BLUE//'INITIALIZING DONE','✔'  
      write (*, *) ANSI_BLUE//'========================================================'//ANSI_RESET
      write (*,'(A,f15.6,A,f15.6,A)') ANSI_BLUE//'Welfare at ss at initialization: V1 = ',v1(pos_zeroA1,pos_zeroA2,pos_zeroOm,pos_zero1,pos_zero2),&
      "; V2 = ",v2(pos_zeroA1,pos_zeroA2,pos_zeroOm,pos_zero1,pos_zero2),ANSI_RESET
      write (*, *) ANSI_BLUE//'========================================================'
      write (*, *) ANSI_BLUE//'SOLVING AUTARKY'
      call value_autarky()
      write (*, *) ANSI_BLUE//'SOLVING AUTARKY DONE ✔'
      write (*, *) ANSI_BLUE//'========================================================'//ANSI_RESET
      write (*,*) ANSI_YELLOW//"Value functions at autarky: ", v1_aut(pos_zeroA1,pos_zeroA2,pos_zero1,pos_zero2), v2_aut(pos_zeroA1,pos_zeroA2,pos_zero1,pos_zero2),ANSI_RESET

      ! call core_model(1.0_wp, 1.0_wp, C1, C2, L1, L2, p1, Q, Y1, res, pred, 1.0_wp,1.0_wp)

!! SOLVE FOR CM ALLOCATION
      ! scale_value=1.0_wp!/abs(v1(pos_zeroOm,pos_zeroA1,pos_zeroA2,pos_zeroK,pos_zero1,pos_zero2))


      write (*, *) ANSI_BLUE//'========================================================'
      write (*, *) ANSI_BLUE//'SOLVING FOR ALLOCATIONS AND OPTIMAL KAPPA'
      write (*, *) ANSI_BLUE//'========================================================'//ANSI_RESET

  
      ! Find optimal kappa using bisection
      call find_zero_arrow_bisection(kappa, min_kappa_res(cnt_mom))
      opt_kappa(cnt_mom) = kappa

      write (*, '(a)') ANSI_YELLOW//'##################################################################'//ANSI_RESET
      write (*, '(a,e20.11)') ANSI_YELLOW//'Found optimal kappa: '//ANSI_RESET, opt_kappa(cnt_mom)
      write (*, '(a,e20.11)') ANSI_YELLOW//'Arrow securities residual: '//ANSI_RESET, min_kappa_res(cnt_mom)
      write (*, '(a)') ANSI_YELLOW//'##################################################################'//ANSI_RESET

      if (.true.) then ! muted for debugging ! BETTER TO COMPUTE THIS EX-POST AFTER BEING SURE OF ACCURACY
         write (*, *) ANSI_BLUE//'========================================================'
         write (*, *) ANSI_BLUE//'COMPUTING WELFARE'
         write (*, *) ANSI_BLUE//'========================================================'//ANSI_RESET

        call toc()


         write (*, *) ANSI_BLUE//'========================================================'
         write (*, *) ANSI_BLUE//'COMPUTING PCU'
         write (*, *) ANSI_BLUE//'========================================================'//ANSI_RESET

         call welfare(welf_mom(cnt_mom, :), opt_kappa(cnt_mom))

         parin0(1:2)=welf_mom(cnt_mom, 1:2)
         if(.false.)then ! which kind of pcu: augmenting stochastic consumption or ss
            ! call find_PCU_multiplier(PCU_s, parin0)
            ! cons_aut_mean(cnt_mom,:)=PCU_s
         else
            pcu_c1_au=C1v_aut(pos_zeroA1,pos_zeroA2,pos_zero1,pos_zero2)!USE SS TO AVOID DEPENDENCE ON MOMENTS !consumption1_autarky()
            ! print*,'pcu_c1',pcu_c1_au

            pcu_c1_au=pcu_CRRA(welf_mom(cnt_mom, 1)-welf_mom(cnt_mom, 3),pcu_c1_au)
            ! print*,'pcu_c1',pcu_c1_au
            ! PRINT*,"---------------------"
            cons_aut_mean(cnt_mom,1)=pcu_c1_au
            pcu_c1_au=C2v_aut(pos_zeroA1,pos_zeroA2,pos_zero1,pos_zero2)
            pcu_c1_au=pcu_CRRA(welf_mom(cnt_mom, 2)-welf_mom(cnt_mom, 4),pcu_c1_au)
            ! print*,'pcu_c1',pcu_c1_au
            ! PRINT*,"---------------------"
            cons_aut_mean(cnt_mom,2)=pcu_c1_au
         endif

         write (*, '(a)') ANSI_BLUE//'##################################################################'//ANSI_RESET
         write (*, '(a)') ANSI_GREEN//'WELFARE of COUNTRY ==>'//cntry//ANSI_RESET
         write (*, '(a)') ANSI_GREEN//'WITH BASE COUNTRY ==>'//names_countries(base_cntry)//ANSI_RESET
         write (*, '(a,5x,a,5x,a,5x,a)') ANSI_BLUE//'Home CM'//ANSI_RESET,ANSI_BLUE//'Foreign CM'//ANSI_RESET,ANSI_BLUE//'Home AU'//ANSI_RESET,ANSI_BLUE//'Foreign AU'//ANSI_RESET
         write (*, '(4e20.8)') welf_mom(cnt_mom, :)
         write (*, '(a)') ANSI_YELLOW//'At kappa=1'
         write (*, '(f15.6,f15.6,a)') v1(pos_zeroA1,pos_zeroA2,pos_zeroOm,pos_zero1,pos_zero2), v2(pos_zeroA1,pos_zeroA2,pos_zeroOm,pos_zero1,pos_zero2),ANSI_RESET
         write (*, '(a)') ANSI_GREEN//'WELFARE GAIN'//ANSI_RESET
         write (*, '(a,5x,a)') ANSI_BLUE//'Home GAIN'//ANSI_RESET, ANSI_BLUE//'Foreign GAIN'//ANSI_RESET
         write (*, '(e20.8,e20.8)') welf_mom(cnt_mom, 1) - welf_mom(cnt_mom, 3), welf_mom(cnt_mom, 2) - welf_mom(cnt_mom, 4)
         write (*, '(a,5x,a,5x)') ANSI_BLUE//'GAIN PCU Home'//ANSI_RESET,ANSI_BLUE//'GAIN PCU Foreign'//ANSI_RESET
         write (*, '(es20.8,es20.8)') cons_aut_mean(cnt_mom,:)

         write (*, '(a)') ANSI_GREEN//'With kappa'//ANSI_RESET
         write (*, '(e20.8)') opt_kappa(cnt_mom)
         write (*, *) ANSI_RED//'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
         write (*, *) ANSI_YELLOW//' Population size (A,B)', n, (1.0-n)
         write (*, *) ANSI_YELLOW//' Domestic Goods Weight (A,B)', nu, nu_b
         write (*, *) ANSI_RED//'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
         write (*, '(a)') ANSI_BLUE//'##################################################################'//ANSI_RESET
      end if

      ! ! SAVE THE SOLUTION TO FILE values and coefficients

      call toc()
      ! direnq = 'mkdir '//trim(dir2)
      ! call execute_command_line('mkdir '//trim(dir))
      call save_results_CRRA(dir, cntry)
   write (*, *) ANSI_BLUE//'========================================================'
   write (*, '(A)') 'Done Country'//trim(cntry)//' ✔'
   call execute_command_line("date")
   write (*, '(A)') '========================================================'//ANSI_RESET
   end do loop_moments! END LOOPING OVER MOMENTS

   header = ['kappa         ', &
   'mean 1        ', 'mean 2        ', &
   'stdev 1       ', 'stdev 2       ', &
   'skewness 1    ', 'skewness 2    ', &
   'kurtosis 1    ', 'kurtosis 2    ', &
   'residual      ', &
   'Welfare 1     ', 'Welfare 2     ', &  ! ← 14 characters
   'Welfare 1 Aut ', 'Welfare 2 Aut ', &
   'Gain PCU 1    ', 'Gain PCU 2    ', &
   'pop_A         ', 'pop_B         ']

   welfout(1:n_moms, 1) = opt_kappa(1:n_moms)
   welfout(1:n_moms, 2) = mom_val_1(1:n_moms, 1)
   welfout(1:n_moms, 3) = mom_val_2(1:n_moms, 1)
   welfout(1:n_moms, 4) = mom_val_1(1:n_moms, 2)
   welfout(1:n_moms, 5) = mom_val_2(1:n_moms, 2)
   welfout(1:n_moms, 6) = mom_val_1(1:n_moms, 3)
   welfout(1:n_moms, 7) = mom_val_2(1:n_moms, 3)
   welfout(1:n_moms, 8) = mom_val_1(1:n_moms, 4)
   welfout(1:n_moms, 9) = mom_val_2(1:n_moms, 4)
   welfout(1:n_moms, 10) = min_kappa_res(1:n_moms)
   welfout(1:n_moms, 11) = welf_mom(1:n_moms, 1)
   welfout(1:n_moms, 12) = welf_mom(1:n_moms, 2)
   welfout(1:n_moms, 13) = welf_mom(1:n_moms, 3)
   welfout(1:n_moms, 14) = welf_mom(1:n_moms, 4)
   welfout(1:n_moms, 15) =cons_aut_mean(1:n_moms,1)
   welfout(1:n_moms, 16) =cons_aut_mean(1:n_moms,2)
   welfout(1:n_moms, 17) = pops(1:n_moms, 1)
   welfout(1:n_moms, 18) = pops(1:n_moms, 2)
   call write_csv(filename=trim(dir)//'//'//trim("kappa_data_"//trim(mom2loop)//".csv"), header=header, data=welfout(1:n_moms,:),&
      rownames=names_countries(pos_subsample(1:n_moms)))
call deallocate_globals()
   write (*, *) ANSI_BLUE//'========================================================'
   write (*, '(A)') 'END OF SIMULATION ✔'
   
call execute_command_line("date")
write (*, '(A)') ANSI_BLUE//'========================================================'//ANSI_RESET
end program
