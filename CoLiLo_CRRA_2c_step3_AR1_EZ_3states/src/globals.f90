!##############################################################################
! MODULE globals
! @ Corsetti-Lipinska-Lombardo
! @ Global parameters/arrays and key functions
!##############################################################################
module globals
   use my_kinds_mod, only: wp
   use bspline_oo_module


   implicit none
!##############################################################################
!                           PARAMETERS OF THE MODEL
!                                   AND
!                      TUNING ONF THE SOLUTION ALGORITHM
!##############################################################################
! SHOULD RECOMPUTE INITIAL CONDITIONS FROM GUESS (TRUE) OR FROM SAVED RESULTS (FALSE)
!##############################################################################
   logical, parameter ::     doinitval = .false.
!##############################################################################
!##############################################################################
   logical, parameter :: verbose = .true. ! whether to print each iteration
!##############################################################################

! ############################################################################
! ################      HOW TO SOLVE THE EULER EQUATIONS ####################
   ! solving (fzero) makes sense only if the value function affects allocations, e.g. EZ. Not with CRRA
   character(len=20), parameter         ::    inner_loop ='update'!'fzero'!'fzero'!  
   logical :: cubic_spline=.true.
!##############################################################################
!##########################  chose whether to use fzero (.true.) or .bisection  #############################
!##############################################################################

!                           ACCURACY OF SOLUTION
!##############################################################################
   real(kind=wp), parameter :: sig = 1.0e-8_wp
   real(kind=wp), parameter :: sig_arrow = 1.0e-10_wp
   real(kind=wp), parameter :: sig_pf = 1.0e-10_wp
   real(kind=wp), parameter :: sig_kappa = sig_pf !1.0e-7_wp
!##############################################################################
!##############################################################################
!                           FACTOR OF ASYMMETRY (1+factor_variance)
!##############################################################################
   real(kind=wp), parameter :: factor_variance = 4.0_wp
!##############################################################################
   ! model parameters
   real(kind=wp), parameter :: gamma = 0.25_wp ! See Rudebush and Swansson on changing sign with 1/gamma>1
   real(kind=wp), parameter :: beta = 0.98_wp
   real(kind=wp), parameter :: lamb = 0.3_wp
   real(kind=wp), parameter :: varphi = 1.75_wp
   real(kind=wp) :: pop1, pop2 ! population level
   real(kind=wp) :: n = 0.5_wp ! this is overwritten in main.f90 based on data
   real(kind=wp):: nu = 0.5_wp! (1.0_wp - (1.0_wp - n)*lamb) !HOME BIAS OF COUNTRY 1 OVERWRITTEN IN main.f90
   real(kind=wp) :: nu_b = 0.5_wp! n*lamb ! ditto NOTE THAT nu_s APPIES TO COUNTRY 1 AND 1-nu_s TO COUTNRY 2
   real(kind=wp), parameter :: price_up = 10.0_wp
   real(kind=wp) :: price_down = 0.005_wp!(1.0e0_wp-nu)**(1/(trade_elast-1.0e0_wp))+(1.0e-8_wp) ! ditto
   real(kind=wp), parameter :: alphha = 0.30_wp
   real(kind=wp), parameter :: chhi = 1.0_wp ! set it to zero if alphha=1
   real(kind=wp), parameter :: delta = 0.019_wp
   real(kind=wp), parameter :: trade_elast = 1.5_wp ! ABOVE 100 THE EQUILIBRIUM SWITCHES TO PERFECT SUBSTITUTION pf=ph
   real(kind=wp), parameter :: theta = 3.5_wp !5.0_wp
   real(kind=wp), parameter :: k0 = 79_wp
   real(kind=wp)            :: rho = 0.8246_wp ! COUNTRY SPECIFIC
   real(kind=wp)            :: rho_A1 = 0.0_wp ! to be reset in main.f90
   real(kind=wp)            :: rho_A2 = 0.0_wp ! to be reset in main.f90
   ! Number of coutries in sample
   integer, parameter :: N_countries = 372
   ! NB: check dimension(XXX) to match subsample elements
   character(len=3), parameter, dimension(2) :: subsample = [ "DEU","MEX"] ![ "DEU","MEX","BRA","THA"], "TUR", "ESP", "KOR", "CHN", "IND", "USA"] ! from the distro of PCU with RoW
   integer, parameter :: actual_cntry = size(subsample,1)
   integer, dimension(actual_cntry) :: pos_subsample
   character(len=2), parameter:: region = '_D' ! possible regions _R=RoW; _E=EME; _A=AE; _D=DEU; _U=USA
   ! variances of stockastic processes
   real(kind=wp):: p1_ar, p2_ar, p3_ar
   real(kind=wp), parameter :: sigma_eps = 0.009499_wp**2 ! minimum variance of the shock from data
   real(kind=wp) :: sigma_eps2
   real(kind=wp) :: sigma_eps3
   real(kind=wp):: m1
   real(kind=wp):: m2
   real(kind=wp) :: m3
   integer,parameter :: deg_int=3
   ! the shock process
   integer, parameter :: NS = 9 ! must be odd to be symmetric around 0
   real(kind=wp) :: pi1(NS, NS), eta1(NS), Cpimx1(NS, NS), Cpimx2(NS, NS)
   real(kind=wp) :: pi2(NS, NS), eta2(NS)
   integer :: pos_zeroA1,pos_zeroA2,pos_zero1, pos_zero2, pos_zeroK,pos_zeroOm! closest to steady state value


! arrays to import markov process
   character(len=5) :: names_countries(N_countries)
   integer, parameter:: n_cols_ar = 11 ! COLUMNS OF THE ar_data_pwt.cvs (MUST BE RIGHT: CHECK THE FILE)
   real(wp), dimension(N_countries, n_cols_ar - 1) :: ar_stats
   real(wp), dimension(N_countries, NS*NS) :: allP
   real(wp), dimension(N_countries, NS):: allX
   integer :: base_cntry !position of base_cntry in name_countries

   ! policy functions FUNCTION OF KAPPA AND TWO EXOGENOUS SHOCKS
   ! integer, parameter :: NK1 = 20 !  must be even (start from 0) to encompas 1
   integer, parameter :: NK1 = 1 !  This is not used anymore, except for my spline interpolation (which I don't use anymore)
   integer, parameter :: NK2 = 6 ! dim of AR terms 
   integer, parameter :: NK3 = 6 ! dim of AR terms 
   integer, parameter :: NK4 = 10 ! dim Omega

   real(wp), allocatable :: tx_v1(:,:,:), ty_v1(:,:,:), tz_v1(:,:,:), tq_v1(:,:,:), bcoef4_v1(:,:,:,:,:)
   real(wp), allocatable :: tx_v2(:,:,:), ty_v2(:,:,:), tz_v2(:,:,:), tq_v2(:,:,:), bcoef4_v2(:,:,:,:,:)

! Policy functions for all variables under optimal kappa
! real(kind=wp) :: kappaS(0:NK1)
   real(kind=wp), allocatable :: C1v_aut(:,:,:,:), C2v_aut(:,:,:,:)
   real(kind=wp), allocatable :: L1v_aut(:,:,:,:), L2v_aut(:,:,:,:)
   real(kind=wp), allocatable :: Y1v_aut(:,:,:,:), Y2v_aut(:,:,:,:), Qv_aut(:,:,:,:)
   real(kind=wp), allocatable :: pfv_aut(:,:,:,:), phv_aut(:,:,:,:)
   real(kind=wp), allocatable :: v1_aut(:,:,:,:), v1_aut_new(:,:,:,:), pf_sols_aut(:,:,:,:)
   real(kind=wp), allocatable :: v2_aut(:,:,:,:), v2_aut_new(:,:,:,:)
   real(kind=wp), allocatable :: s1(:,:,:,:,:), s1_new(:,:,:,:,:)
   real(kind=wp), allocatable :: A1(:), A2(:), OmegaS(:)

   real(kind=wp), parameter :: Omega_l = 0.90_wp ! Difference between MUC
   real(kind=wp), parameter :: Omega_u = 1.10_wp ! must consider the resource constraint: C1<(n*A1+(1-n)*A2)/n
   real(kind=wp), parameter :: A1_l = 0.60_wp, A2_l = A1_l, A1_u = 1.40_wp, A2_u = A1_u ! these will be overwritten in main.f90 based on data
   real(kind=wp),  dimension(3) ::    LB_all = [A1_l, A2_l,Omega_l], UB_all = [A1_u, A2_u,Omega_u] !must be updated later
    
   ! variables to numerically determine policy function

   real(kind=wp) :: con_lev_v1, con_lev_v2, con_lev_s1, x_in1, x_in2, x_in3
   real(kind=wp), dimension(2) :: x_input ! three states and kappa
   logical :: check
! COLORS
   character(len=*), parameter :: ANSI_RESET = achar(27)//"[0m"
   character(len=*), parameter :: ANSI_BACKGROUND = achar(27)//"[102m" !light green
   character(len=*), parameter :: ANSI_RED = achar(27)//"[31m"
   character(len=*), parameter :: ANSI_GREEN = achar(27)//"[32m"
   character(len=*), parameter :: ANSI_BLUE = achar(27)//"[1;36m" !cycan which reads on black; blue is too dark
   character(len=*), parameter :: ANSI_YELLOW = achar(27)//"[1;93m"
   ! variables to communicate with function
   real(kind=wp) :: k_com, parin(7)
   integer :: is_com
   ! ##### PARAMETER VECTOR FOR OPTIMIZATION PARALLELIZED
   ! real(kind=wp) :: parin0All(7, 500) !# assuming have no more than 500 processors

   

   ! computation of the Jacobian
   !##############################################################################
   !# looping on kappas
   !#################################################################################
   integer, parameter :: nkappas = 1000
   real(kind=wp), parameter :: step_kappa = 0.1_wp, kappa_min = 0.8_wp, kappa_max = 1.2_wp
   real(kind=wp), dimension(nkappas + 1) :: possible_kappas
   real(kind=wp), dimension(nkappas) :: resid_A0
   ! matrices to save results: values and coefficients
   real(kind=wp), allocatable ::v1(:,:,:,:,:), v2(:,:,:,:,:), v1_test(:,:,:,:,:), v1_new(:,:,:,:,:), v2_new(:,:,:,:,:), pf_sols(:,:,:,:,:), phv(:,:,:,:,:)
   real(kind=wp), allocatable ::Ev1(:,:,:,:,:), Ev2(:,:,:,:,:)
   real(kind=wp), allocatable :: Ev1m1(:,:,:,:,:), Ev2m1(:,:,:,:,:)
   real(kind=wp), allocatable :: C1v(:,:,:,:,:), C2v(:,:,:,:,:), Y1v(:,:,:,:,:), L1v(:,:,:,:,:), L2v(:,:,:,:,:)
   real(kind=wp), allocatable :: Es1(:,:,:,:,:), OmegaSplusv(:,:,:,:,:)
   real(kind=wp), allocatable :: coeff_v1(:,:,:,:,:), coeff_v2(:,:,:,:,:), coeff_s1(:,:,:,:,:)
   real(kind=wp), allocatable :: matrix_csv_v1(:,:), matrix_csv_v2(:,:), matrix_csv_s1(:,:), matrix_csv_pf(:,:)
   real(kind=wp), allocatable :: matrix_csv_aut_v1(:,:), matrix_csv_aut_v2(:,:), matrix_csv_aut_pf(:,:)
   real(kind=wp), allocatable :: matrix_coefs_csv_v1(:,:), matrix_coefs_csv_v2(:,:), matrix_coefs_csv_s1(:,:)
   real(kind=wp) :: meanval_v1, meanval_v2
   real(kind=wp), allocatable :: resh_dummy(:)

   !##############################################################################
   !# END looping on kappas
   !#################################################################################

   !##############################################################################
   !# looping on moments parameters
   !#################################################################################
   INTEGER, parameter :: n_moms_ = 4000 ! this is max number of runs
   Integer :: n_moms ! this is the specific number of runs
   integer, dimension(n_moms_)                     :: opt_kappa_pos
   real(kind=wp), dimension(n_moms_)   :: opt_kappa, min_kappa_res
   real(kind=wp), dimension(n_moms_, 4) :: mom_val_1, mom_val_2
   real(kind=wp), dimension(n_moms_, 8) :: mom_range
   real(kind=wp), dimension(n_moms_) ::  ToT_moms,ToT_aut_moms,pkh_moms,pkf_moms,pkh_aut_moms,pkf_aut_moms,Qv_aut_moms, Qv_moms
   real(kind=wp), dimension(n_moms_,2) :: cons_aut_mean
   real(kind=wp)                      :: mom_val_addition
   real(kind=wp)  :: inputpar1(9), inputpar2(9)!! parameter for discretization
   character*1 :: mom2loop
   !##############################################################################
   !#  end looping on moments parameters
   !#################################################################################

   ! initialize random number generators
   logical, parameter :: tbox_seed = .true., fixed = .true.

   ! numerical parameter

   integer, parameter :: itermax = 10000
   real(kind=wp), parameter :: smooth = 1.0_wp, smooth_outer_s = smooth, smooth_pf = smooth ! smoothing parameter of FPI update
   real(kind=wp), parameter :: smooth_alittle = 0.1_wp, smooth_alot = 0.1_wp
   
   real(kind=wp) :: smooth_outer = 0.3_wp
   real(kind=wp),parameter :: smooth_arrow = 1.0_wp,smooth_aut = 0.3_wp

   ! counter variables
   integer :: it, iter

   ! time path of consumption and capital
   integer, parameter :: TT = 5000
   real(kind=wp) :: c1_t(0:TT), c2_t(0:TT),L1_t(0:TT), L2_t(0:TT), eta1_t(0:TT), eta2_t(0:TT), v1_t(0:TT), v2_t(0:TT), Omega_t(0:TT), mu_t(0:TT)
   real(kind=wp) ::c1_aut_t(0:TT), c2_aut_t(0:TT), L1_aut_t(0:TT), L2_aut_t(0:TT), v1_aut_t(0:TT), v2_aut_t(0:TT)
   !##############################################################################
   !# simulation for looping on kappas
   !#################################################################################
   integer :: is_t(0:TT), is1_t(0:TT), is2_t(0:TT)
   integer, dimension(TT) :: totindex1, totindex2
   real(kind=wp), dimension(TT) :: eta_simul1, eta_simul2
   real(kind=wp), dimension(TT) :: zindex1, zindex2
   !##############################################################################
   !# END simulation for looping on kappas
   !#################################################################################

   real(kind=wp), allocatable :: coeff_v1_aut(:,:,:,:), coeff_v2_aut(:,:,:,:)

   interface
        subroutine core_model(pf, kappa, OmegaS_in, C, Cs, L, Ls, ph, Q, Y1, res, pred, D, Ds)
         use my_kinds_mod, only: wp
         implicit none
            real(kind=wp), intent(in):: kappa, pf, D, Ds, OmegaS_in
         real(kind=wp)            :: mu, mus, YW
         integer           :: cntr
         real(kind=wp), intent(out) :: C, Cs, L, Ls, ph, Q, pred, res, Y1
      end subroutine core_model
   end interface

   ! interface
   ! function kappafunc(x_input, parin0)
   !     use my_toolbox, only: spline_eval
   !     use my_kinds_mod, only: wp
   !     implicit none
   !     real(kind=wp), intent(in) :: x_input
   !     real(kind=wp) :: kappafunc
   !     real(kind=wp), intent(in) :: parin0(:)

   ! end function
   ! end interface
   interface
      subroutine core_model_autarky(pf, C, Cs, L, Ls, ph, Q, Y1, res, pred, D, Ds)
         use my_kinds_mod, only: wp
         implicit none
         real(kind=wp):: pf
         real(kind=wp), intent(in)             ::  D, Ds
         real(kind=wp), intent(out) :: C, Cs, L, Ls, ph, Q, pred, res, Y1
      end subroutine core_model_autarky
   end interface
   interface
      subroutine write_csv_portable(filename, data)
         use my_kinds_mod
         character(len=*), intent(in) :: filename
         real(kind=wp), intent(in)    :: data(:, :)
      end subroutine write_csv_portable
   end interface
   ! interface
   !     pure function pcu_CRRA(W, EC, EL) result(pcu)
   !         use my_kinds_mod
   !         implicit none

   !         real(kind=wp), intent(in) :: W, EC, EL
   !         real(kind=wp) :: pcu
   !     end function
   ! end interface

   interface
      function asset_prices(kappa) result(ap)
         use my_kinds_mod
         implicit none

         real(kind=wp) :: C1, C2, Q, L1, L2, ph, pf
         real(kind=wp) :: Y1, Y2, ap(2)
         real(kind=wp) :: Eap1, Eap2
         real(kind=wp), intent(in)::  kappa
         integer :: is1_n, is2_n,  is_p1, is_p2
      end function
   end interface

   interface
      function asset_prices_aut() result(ap)
         use my_kinds_mod
         implicit none

         real(kind=wp) :: C1, C2, Q, L1, L2, ph, pf
         real(kind=wp) :: Y1, Y2, ap(2)
         real(kind=wp) :: Eap1, Eap2
         integer :: is1_n, is2_n,  is_p1, is_p2
      end function
   end interface

   ! interface
   ! function solve_pf_autarky(x_in, parin01) result(result_)

   !     implicit none
   ! ! must use kind 8 to be compatible with funcv in toolbox
   !     real(wp), intent(in) :: x_in
   !     real(wp), intent(in) :: parin01(:)
   !     real(wp) :: result_,D,Ds,C,Cs,L,Ls,ph,Q,Y1,pred
   !     end function solve_pf_autarky
   ! end interface
contains
   subroutine solve_value(nres, size_x, x_input, grad, need_gradient, parin0)
      implicit none
      integer, intent(in) :: size_x
      real(kind=wp), intent(inout) :: x_input(size_x)
      real(kind=wp) :: res(size_x), temp(size_x),nres, grad(size_x)
      real(kind=wp), intent(in):: parin0(100)
      integer::need_gradient
      temp = x_input


      call myfunc_update(size_x, x_input, res, parin0)

      nres=norm2(res,1)
   end subroutine
   ! the first order condition
   function foc2(x_input, parin0)
      implicit none
      real(kind=wp), intent(in) :: x_input(:)
      real(kind=wp) :: foc2(size(x_input, 1)), temp(size(x_input, 1))
      real(kind=wp), intent(in) :: parin0(:)
      temp = x_input
      call myfunc_update(size(x_input, 1), temp, foc2, parin0)
   end function



   ! solver for kappa
   function kappafunc(xkappa, parin0)
      use my_toolbox, only: spline_eval
      use my_kinds_mod, only: wp
      implicit none
      real(kind=wp), intent(in) :: xkappa
      real(kind=wp) :: kappafunc
      real(kind=wp), intent(in) :: parin0(:)
      kappafunc = spline_eval([A1(pos_zeroA1),A2(pos_zeroA2),OmegaS(pos_zeroOm)], coeff_s1(:,:,:, pos_zero1, pos_zero2), LB_all, UB_all)
   end function

   subroutine kappafunc_subr(nres, size_x, xkappa, grad, need_gradient, parin0)
      use my_toolbox, only: spline_eval
      use my_kinds_mod, only: wp
      implicit none
      integer, intent(in) :: size_x
      real(kind=wp), intent(inout) :: xkappa(size_x)
      real(kind=wp) :: results(size(xkappa, 1)), grad(size_x), nres
      real(kind=wp), intent(in) :: parin0(100)
      integer:: need_gradient
      results = spline_eval([A1(pos_zeroA1),A2(pos_zeroA2),OmegaS(pos_zeroOm)], coeff_s1(:,:,:, pos_zero1, pos_zero2), LB_all, UB_all)

      nres = abs(results(1))
   end subroutine
! constraint for prices
   subroutine myconstraint(nres, size_x, x_in, grad, need_gradient, parin0)
      use my_kinds_mod
      implicit none
      integer, intent(in) :: size_x
      real(kind=wp), intent(in) :: x_in(size_x)
      real(kind=wp), intent(in) :: parin0(100)
      real(kind=wp)             ::  pf
      real(kind=wp) ::  grad(size_x), nres
      integer ::  need_gradient
      pf= x_in(1)

      if (need_gradient .ne. 0) then
         grad(1) =-(((1.0e0_wp-nu)*pf**(-trade_elast)*((1.0e0_wp-(1.0e0_wp-nu)*pf**(1.0e0_wp-trade_elast))/nu)**(-1.0e0_wp+1.0e0_wp/(1.0e0_wp-trade_elast)))/nu)
      end if

      nres= ((1.0e0_wp - (1.0e0_wp - nu)*pf**(1.0e0_wp - trade_elast))/nu)**(1.0e0_wp/(1.0e0_wp - trade_elast))
   end subroutine myconstraint

    function func_omega(xin_V1,xin_V2,OmegaS_n,A1_n,A2_n,parin0) result(OmegaSplus)
        use my_kinds_mod
        use my_toolbox,only: spline_eval 
        implicit none
        real(kind=wp), intent(in) :: xin_V1,xin_V2,OmegaS_n,A1_n,A2_n
        real(kind=wp), intent(in) :: parin0(:)
        real(kind=wp) :: v1plusm1,v2plusm1,x_states(3),Omega1,Omega2,OmegaSplus,kappa
        integer :: is_p1,is_p2,ik2_n,ik3_n,ik4_n, is1_n,is2_n
        is1_n = INT(parin0(1)) ! current shock index
        is2_n = INT(parin0(2)) ! current shock index
        kappa = parin0(3)
        ik2_n = INT(parin0(4))
        ik3_n = INT(parin0(5))
        ik4_n = INT(parin0(6))
            ! compute E_{t-1}V_t, which depends on state t-1 and expected innovations    
        Ev1m1(ik2_n, ik3_n, ik4_n, is1_n,is2_n) = 0e0_wp
        Ev2m1(ik2_n, ik3_n, ik4_n, is1_n,is2_n) = 0e0_wp
        x_states = [min(A1_n, UB_all(1)), min(A2_n, UB_all(2)), min(OmegaS_n, UB_all(3))]
        do is_p1 = 1, NS
            do is_p2 = 1, NS
                ! lagged expectations
                v1plusm1 = spline_eval(x_states, coeff_v1(:,:,:, is_p1, is_p2), LB_all, UB_all)
                v2plusm1 = spline_eval(x_states, coeff_v2(:,:,:, is_p1, is_p2), LB_all, UB_all)
                if(v1plusm1<0)print*,ANSI_RED//"Warning: negative value in func_omega v1plusm1"

                Ev1m1(ik2_n, ik3_n, ik4_n,is1_n,is2_n) = Ev1m1(ik2_n, ik3_n, ik4_n, is1_n,is2_n) &
                + pi1(is1_n, is_p1)*pi2(is2_n, is_p2)*(v1plusm1)**theta ! See Rudebush and Swansson on changing sign with 1/gamma>1
                Ev2m1(ik2_n, ik3_n, ik4_n, is1_n,is2_n) = Ev2m1(ik2_n, ik3_n, ik4_n, is1_n,is2_n) &
                + pi1(is1_n, is_p1)*pi2(is2_n, is_p2)*(v2plusm1)**theta

            end do
        end do

        ! compute allocations at time t
        Omega1 = beta*(Ev1m1(ik2_n, ik3_n, ik4_n, is1_n,is2_n))**(1.0_wp/theta - 1.0_wp)*(-xin_V1)**(theta - 1.0_wp) ! See RUdebusch and Swansson for change of sign
        Omega2 = beta*(Ev2m1(ik2_n, ik3_n, ik4_n, is1_n,is2_n))**(1.0_wp/theta - 1.0_wp)*(-xin_V2)**(theta - 1.0_wp)

        ! update state variables
        ! OmegaS (ie get OmegaS_t, which enters consumption shares and thus defines allocations)
        OmegaSplus = Omega2/Omega1*OmegaS_n
    end function
    subroutine solve_pf(nres, size_x, x_in, grad, need_gradient, parin01)
        use my_kinds_mod
        implicit none
        integer, intent(in) :: size_x
        REAL(kind=wp), INTENT(INOUT), DIMENSION(size_x) :: x_in
        real(kind=wp), INTENT(IN), DIMENSION(100) :: parin01
        real(kind=wp)             ::   C1, C2, L1, L2,  ph, Q, Y1, Y2,xin_V1,xin_V2,OmegaS_n,A1_n,A2_n,OmegaSplus,A1plus,A2plus
        real(kind=wp) :: res(size_x), pred, kappa, grad(size_x), nres,pf
        integer :: is1_n, is2_n, ik1_n,ik2_n,ik3_n, ik4_n, need_gradient

        is1_n = INT(parin01(1)) ! current shock index
        is2_n = INT(parin01(2)) ! current shock index
        kappa = (parin01(3))
        ik2_n = INT(parin01(4))
        ik3_n = INT(parin01(5))
        ik4_n = INT(parin01(6))

        xin_V1=parin01(8)
        xin_V2=parin01(9)
        A1_n=parin01(10)
        A2_n=parin01(11)
        OmegaS_n = OmegaS(ik4_n)



        OmegaSplus=func_omega(xin_V1,xin_V2,OmegaS_n,A1_n,A2_n,parin01(1:7))
! OmegaSplus=1.0_wp ! for debugging

        A1plus = A1_n**rho_A1*exp(eta1(is1_n))
        A2plus = A2_n**rho_A2*exp(eta2(is2_n))

        pf=x_in(1)
         call core_model(pf=pf, kappa=kappa, OmegaS_in=OmegaSplus, C=C1, Cs=C2, L=L1, Ls=L2, &
         ph=ph, Q=Q, Y1=Y1, res=nres, pred=pred, D=A1plus, Ds=A2plus)
        
    end subroutine



    subroutine solve_pf_minpack(size_x, x_in, res,iflag, parin01)
        use my_kinds_mod
        implicit none
        integer, intent(in) :: size_x
        real(kind=wp), intent(in) :: x_in(size_x)
        real(kind=wp),intent(in),optional,dimension(:) :: parin01
        integer,intent(inout) :: iflag
        real(kind=wp) , intent(out):: res(size_x)
        real(kind=wp)             ::   C1, C2, L1, L2,  ph, Q, Y1, Y2,xin_V1,xin_V2,OmegaS_n,A1_n,A2_n,OmegaSplus,A1plus,A2plus
        real(kind=wp) :: pred, kappa, grad(size_x), nres,pf
        integer :: is1_n, is2_n, ik2_n,ik3_n, ik4_n, need_gradient

        is1_n = INT(parin01(1)) ! current shock index
        is2_n = INT(parin01(2)) ! current shock index
        kappa = parin01(3)
        ik2_n = INT(parin01(4))
        ik3_n = INT(parin01(5))
        ik4_n = INT(parin01(6))

        xin_V1=parin01(8)
        xin_V2=parin01(9)
        A1_n=parin01(10)
        A2_n=parin01(11)
        
        OmegaS_n = OmegaS(ik4_n)

    

        OmegaSplus=func_omega(xin_V1,xin_V2,OmegaS_n,A1_n,A2_n,parin01(1:7))
        ! OmegaSplus=1.0_wp ! for debugging
        A1plus = A1_n**rho_A1*exp(eta1(is1_n))
        A2plus = A2_n**rho_A2*exp(eta2(is2_n))
        pf=x_in(1)

         call core_model(pf=pf, kappa=kappa, OmegaS_in=OmegaSplus, C=C1, Cs=C2, L=L1, Ls=L2, &
         ph=ph, Q=Q, Y1=Y1, res=nres, pred=pred, D=A1plus, Ds=A2plus)
        res(1)=nres
    end subroutine
    subroutine allocate_globals()
      
      implicit none

      allocate(v1(0:NK2,0:NK3,0:NK4,1:NS,1:NS))
      allocate(v2(0:NK2,0:NK3,0:NK4,1:NS,1:NS))
      allocate(v1_test(0:NK2,0:NK3,0:NK4,1:NS,1:NS))
      allocate(v1_new(0:NK2,0:NK3,0:NK4,1:NS,1:NS))
      allocate(v2_new(0:NK2,0:NK3,0:NK4,1:NS,1:NS))
      allocate(pf_sols(0:NK2,0:NK3,0:NK4,1:NS,1:NS))
      allocate(phv(0:NK2,0:NK3,0:NK4,1:NS,1:NS))
      allocate(Ev1(0:NK2,0:NK3,0:NK4,1:NS,1:NS))
      allocate(Ev2(0:NK2,0:NK3,0:NK4,1:NS,1:NS))
      allocate(Ev1m1(0:NK2,0:NK3,0:NK4,1:NS,1:NS))
      allocate(Ev2m1(0:NK2,0:NK3,0:NK4,1:NS,1:NS))
      allocate(C1v(0:NK2,0:NK3,0:NK4,1:NS,1:NS))
      allocate(C2v(0:NK2,0:NK3,0:NK4,1:NS,1:NS))
      allocate(Y1v(0:NK2,0:NK3,0:NK4,1:NS,1:NS))
      allocate(L1v(0:NK2,0:NK3,0:NK4,1:NS,1:NS))
      allocate(L2v(0:NK2,0:NK3,0:NK4,1:NS,1:NS))
      allocate(Es1(0:NK2,0:NK3,0:NK4,1:NS,1:NS))
      allocate(OmegaSplusv(0:NK2,0:NK3,0:NK4,1:NS,1:NS))

      allocate(matrix_csv_v1((NK2 + 1)*(NK3 + 1)*(NK4 + 1)*NS, NS))
      allocate(matrix_csv_v2((NK2 + 1)*(NK3 + 1)*(NK4 + 1)*NS, NS))
      allocate(matrix_csv_s1((NK2 + 1)*(NK3 + 1)*(NK4 + 1)*NS, NS))
      allocate(matrix_csv_pf((NK2 + 1)*(NK3 + 1)*(NK4 + 1)*NS, NS))

      allocate(matrix_csv_aut_v1((NK2 + 1)*(NK3 + 1)*NS, NS))
      allocate(matrix_csv_aut_v2((NK2 + 1)*(NK3 + 1)*NS, NS))
      allocate(matrix_csv_aut_pf((NK2 + 1)*(NK3 + 1)*NS, NS))

      allocate(matrix_coefs_csv_v1((NK2 + 3)*(NK3 + 3)*(NK4 + 1)*NS, NS))
      allocate(matrix_coefs_csv_v2((NK2 + 3)*(NK3 + 3)*(NK4 + 1)*NS, NS))
      allocate(matrix_coefs_csv_s1((NK2 + 3)*(NK3 + 3)*(NK4 + 1)*NS, NS))

      allocate(resh_dummy((NK2+1)*(NK3+1)*(NK4 + 1)*NS*NS))

      allocate(coeff_v1(1:NK2 + 3,1:NK3 + 3,1:NK4+3, NS, NS))
      allocate(coeff_v2(1:NK2 + 3,1:NK3 + 3,1:NK4+3, NS, NS))
      allocate(coeff_s1(1:NK2 + 3,1:NK3 + 3,1:NK4+3, NS, NS))

      allocate(tx_v1(NK2+1+deg_int,NS,NS))
      allocate(ty_v1(NK3+1+deg_int,NS,NS))
      allocate(tz_v1(NK4+1+deg_int,NS,NS))
      allocate(bcoef4_v1(0:NK2,0:NK3,0:NK4,NS,NS))
      allocate(tx_v2(NK2+1+deg_int,NS,NS))
      allocate(ty_v2(NK3+1+deg_int,NS,NS))
      allocate(tz_v2(NK4+1+deg_int,NS,NS))
      allocate(bcoef4_v2(0:NK2,0:NK3,0:NK4,NS,NS))
      allocate(C1v_aut(0:NK2,0:NK3,NS,NS))
      allocate(C2v_aut(0:NK2,0:NK3,NS,NS))
      allocate(L1v_aut(0:NK2,0:NK3,NS,NS))
      allocate(L2v_aut(0:NK2,0:NK3,NS,NS))
      allocate(Y1v_aut(0:NK2,0:NK3,NS,NS))
      allocate(Y2v_aut(0:NK2,0:NK3,NS,NS))
      allocate(Qv_aut(0:NK2,0:NK3,NS,NS))
      allocate(pfv_aut(0:NK2,0:NK3,NS,NS))
      allocate(phv_aut(0:NK2,0:NK3,NS,NS))
      allocate(v1_aut(0:NK2,0:NK3,1:NS,1:NS))
      allocate(v1_aut_new(0:NK2,0:NK3,1:NS,1:NS))
      allocate(pf_sols_aut(0:NK2,0:NK3,1:NS,1:NS))
      allocate(v2_aut(0:NK2,0:NK3,1:NS,1:NS))
      allocate(v2_aut_new(0:NK2,0:NK3,1:NS,1:NS))
      allocate(s1(0:NK2,0:NK3,0:NK4,1:NS,1:NS))
      allocate(s1_new(0:NK2,0:NK3,0:NK4,1:NS,1:NS))
      allocate(A1(0:NK2))
      allocate(A2(0:NK3))
      allocate(OmegaS(0:NK4))
      allocate(coeff_v1_aut(1:NK2+3, 1:NK3+3, NS, NS))
      allocate(coeff_v2_aut(1:NK2+3, 1:NK3+3, NS, NS))
   end subroutine allocate_globals

   subroutine deallocate_globals()
      implicit none
      if (allocated(tx_v1)) deallocate(tx_v1)
      if (allocated(ty_v1)) deallocate(ty_v1)
      if (allocated(tz_v1)) deallocate(tz_v1)
      if (allocated(tq_v1)) deallocate(tq_v1)
      if (allocated(bcoef4_v1)) deallocate(bcoef4_v1)
      if (allocated(tx_v2)) deallocate(tx_v2)
      if (allocated(ty_v2)) deallocate(ty_v2)
      if (allocated(tz_v2)) deallocate(tz_v2)
      if (allocated(tq_v2)) deallocate(tq_v2)
      if (allocated(bcoef4_v2)) deallocate(bcoef4_v2)
      if (allocated(C1v_aut)) deallocate(C1v_aut)
      if (allocated(C2v_aut)) deallocate(C2v_aut)
      if (allocated(L1v_aut)) deallocate(L1v_aut)
      if (allocated(L2v_aut)) deallocate(L2v_aut)
      if (allocated(Y1v_aut)) deallocate(Y1v_aut)
      if (allocated(Y2v_aut)) deallocate(Y2v_aut)
      if (allocated(Qv_aut)) deallocate(Qv_aut)
      if (allocated(pfv_aut)) deallocate(pfv_aut)
      if (allocated(phv_aut)) deallocate(phv_aut)
      if (allocated(v1_aut)) deallocate(v1_aut)
      if (allocated(v1_aut_new)) deallocate(v1_aut_new)
      if (allocated(pf_sols_aut)) deallocate(pf_sols_aut)
      if (allocated(v2_aut)) deallocate(v2_aut)
      if (allocated(v2_aut_new)) deallocate(v2_aut_new)
      if (allocated(s1)) deallocate(s1)
      if (allocated(s1_new)) deallocate(s1_new)
      if (allocated(A1)) deallocate(A1)
      if (allocated(A2)) deallocate(A2)
      if (allocated(OmegaS)) deallocate(OmegaS)
      if (allocated(v1)) deallocate(v1)
      if (allocated(v2)) deallocate(v2)
      if (allocated(v1_test)) deallocate(v1_test)
      if (allocated(v1_new)) deallocate(v1_new)
      if (allocated(v2_new)) deallocate(v2_new)
      if (allocated(pf_sols)) deallocate(pf_sols)
      if (allocated(phv)) deallocate(phv)
      if (allocated(Ev1)) deallocate(Ev1)
      if (allocated(Ev2)) deallocate(Ev2)
      if (allocated(Ev1m1)) deallocate(Ev1m1)
      if (allocated(Ev2m1)) deallocate(Ev2m1)
      if (allocated(C1v)) deallocate(C1v)
      if (allocated(C2v)) deallocate(C2v)
      if (allocated(Y1v)) deallocate(Y1v)
      if (allocated(L1v)) deallocate(L1v)
      if (allocated(L2v)) deallocate(L2v)
      if (allocated(Es1)) deallocate(Es1)
      if (allocated(OmegaSplusv)) deallocate(OmegaSplusv)
      if (allocated(matrix_csv_v1)) deallocate(matrix_csv_v1)
      if (allocated(matrix_csv_v2)) deallocate(matrix_csv_v2)
      if (allocated(matrix_csv_s1)) deallocate(matrix_csv_s1)
      if (allocated(matrix_csv_pf)) deallocate(matrix_csv_pf)
      if (allocated(matrix_csv_aut_v1)) deallocate(matrix_csv_aut_v1)
      if (allocated(matrix_csv_aut_v2)) deallocate(matrix_csv_aut_v2)
      if (allocated(matrix_csv_aut_pf)) deallocate(matrix_csv_aut_pf)
      if (allocated(matrix_coefs_csv_v1)) deallocate(matrix_coefs_csv_v1)
      if (allocated(matrix_coefs_csv_v2)) deallocate(matrix_coefs_csv_v2)
      if (allocated(matrix_coefs_csv_s1)) deallocate(matrix_coefs_csv_s1)
      if (allocated(coeff_v1)) deallocate(coeff_v1)
      if (allocated(coeff_v2)) deallocate(coeff_v2)
      if (allocated(coeff_s1)) deallocate(coeff_s1)
      if (allocated(coeff_v1_aut)) deallocate(coeff_v1_aut)
      if (allocated(coeff_v2_aut)) deallocate(coeff_v2_aut)
   end subroutine deallocate_globals
end module
