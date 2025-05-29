module optimization_interface_mod

    use omp_lib
    use :: toolbox, only:fzero

    ! use iso_fortran_env, only: nwrite => output_unit
    ! use ogpf

    implicit none
contains

    subroutine optimize(xin, func_, lb, ub, fdata, code_, constr)
        use my_kinds_mod, only: wp
        use globals
        implicit none
        real(kind=wp), intent(inout) :: xin(:)
        integer :: size_x
        integer, optional, intent(in) :: code_
        integer, optional, intent(in) :: constr
        real(kind=wp) :: lb(size(xin, 1)), ub(size(xin, 1)), tolc(size(xin, 1)), x(size(xin, 1)), tol(size(xin, 1))
        real(kind=wp), parameter::tolly = 1.0e-9_wp
        real(kind=wp) :: fdata(:)
        integer :: solver, slv
        integer*8 :: opt
        real(kind=wp) :: minf
        integer :: ires, nevals
        
      
      

        include 'nlopt.f'

        ! interface for the function
        interface
            subroutine func_(nres, size_x, x_in, grad, need_gradient, parin0)
                use my_kinds_mod
                implicit none
                integer, intent(in) :: size_x
                real(kind=wp), intent(in) :: x_in(size_x)
                real(kind=wp), intent(in) :: parin0(7)
                real(kind=wp)             ::  D, Ds, Dc, C, Cs, Cc, L, Ls, Lc, pA, Qb, Qc, Ya, Yb
                real(kind=wp) :: res(size_x), pred(size_x), kappa(size_x), grad(size_x), nres
                integer :: is1_n, is2_n, is3_n, ik1_n, ik2_n, need_gradient
            end subroutine func_
        end interface

        ! interface
        !     subroutine myconstraint(nres, size_x, x_in, grad, need_gradient, parin0)
        !         use my_kinds_mod
        !         implicit none
        !         integer, intent(in) :: size_x
        !         real(kind=wp), intent(in) :: x_in(size_x)
        !         real(kind=wp), intent(in) :: parin0(7)
        !         real(kind=wp)             ::  pB, pC
        !         real(kind=wp) ::  grad(size_x), nres
        !         integer ::  need_gradient
        !     end subroutine myconstraint
        ! end interface

        x = xin
        size_x = size(xin, 1)

        if (present(code_) .neqv. .true.) then
            write (*, *) "--------------------------"
            write (*, *) "\033[95m CHOOSE SOLVER \033[0m"
            write (*, *) "--------------------------"
            write (*, *) "1=\033[95m NLOPT_LN_NELDERMEAD \033[0m"
            write (*, *) "2=\033[95m NLOPT_GN_CRS2_LM \033[0m"
            write (*, *) "3=\033[95m NLOPT_LN_SBPLX \033[0m"
            write (*, *) "4=\033[95m NLOPT_LD_MMA \033[0m"
       write (*, *) "5=\033[95m NLOPT_GN_DIRECT_L_RAND: DIRECT is the DIviding RECTangles algorithm for global optimization \033[0m"

            write (*, *) "6=\033[95m NLOPT_GN_ESCH \033[0m"
            write (*, *) "7=\033[95m NLOPT_LN_COBYLA \033[0m"
            write (*, *) "8=\033[95m SLSQP: NLOPT_LD_SLSQP \033[0m"
            write (*, *) "9=\033[95m AUGLAG: NLOPT_AUGLAG  \033[0m"
            write (*, *) "10=\033[95m NLOPT_LD_AUGLAG \033[0m"
            write (*, *) "11=\033[95m ISRES: NLOPT_GN_ISRES \033[0m"

            write (*, *) "12=\033[95m  LD_LBFGS: Limited-memory BFGS\033[0m"
            write (*, *) "13=\033[95m LD_TNEWTON_PRECOND_RESTART: Preconditioned truncated Newton with restarts\033[0m"
            write (*, *) "14=\033[95m  LD_TNEWTON_PRECOND: Preconditioned truncated Newton\033[0m"
            write (*, *) "15=\033[95m LD_TNEWTON_RESTART: Truncated Newton with restarts\033[0m"
            write (*, *) "16=\033[95m  LD_TNEWTON: Truncated Newton\033[0m"
            write (*, *) "17=\033[95m  LD_VAR2: Shifted limited-memory variable-metric\033[0m"
            write (*, *) "18=\033[95m  LD_VAR1: Shifted limited-memory variable-metric\033[0m"

! Constrained Optimization:

            write (*, *) "19=\033[95m LD_SLSQP: Sequential Least-Squares Quadratic Programming (SLSQP)\033[0m"

! Global Optimization:

            write (*, *) "20=\033[95m GD_STOGO: Gradient-based global optimization\033[0m"
            write (*, *) "21=\033[95m GD_STOGO_RAND: Randomized gradient-based global optimization\033[0m"
            write (*, *) "22=\033[95m  NLOPT_GN_DIRECT: DIRECT is the DIviding RECTangles algorithm for global optimization \033[0m"

           write (*, *) "23=\033[95m NLOPT_GN_DIRECT_L: DIRECT is the DIviding RECTangles algorithm for global optimization \033[0m"
            write (*, *) "24=\033[95m  NLOPT_GN_DIRECT_L_RAND_NOSCAL: DIRECT is the DIviding RECTangles algorithm for global optimization \033[0m"
      write (*, *) "25=\033[95m NLOPT_GN_DIRECT_NOSCAL: DIRECT is the DIviding RECTangles algorithm for global optimization \033[0m"

    write (*, *) "26=\033[95m NLOPT_GN_DIRECT_L_NOSCAL: DIRECT is the DIviding RECTangles algorithm for global optimization \033[0m"
    write (*,*) "27=\033[95m NLOPT_GN_ORIG_DIRECT: As above but with inequality constraints\033[0m"
    write (*,*) "28=\033[95m NLOPT_GN_ORIG_DIRECT: As above but with inequality constraints\033[0m"
    write (*,*) "29=\033[95m AGS can handle arbitrary objectives and nonlinear inequality constraints. Also bound constraints are required for this method\033[0m"
            write (*, *) "<0 ==>\033[95m  NOT RUNNING SOLVERS \033[0m"
            write (*, *) "--------------------------"
            write (*, *) "--------------------------"
            read (*, '(i)') solver
        else
            solver = code_
        end if
        select case (solver)
        case (1)
            slv = NLOPT_LN_NELDERMEAD
        case (2)
            slv = NLOPT_GN_CRS2_LM
        case (3)
            slv = NLOPT_LN_SBPLX
        case (4)
            slv = NLOPT_LD_MMA
        case (5)
            slv = NLOPT_GN_DIRECT_L_RAND
        case (6)
            slv = NLOPT_GN_ESCH
        case (7)
            slv = NLOPT_LN_COBYLA
        case (8)
            slv = NLOPT_LD_SLSQP
        case (9)
            slv = NLOPT_AUGLAG
        case (10)
            slv = NLOPT_LD_AUGLAG
        case (11)
            slv = NLOPT_GN_ISRES
        case (12)
            slv = NLOPT_LD_LBFGS
            ! LD_LBFGS: Limited-memory BFGS
        case (13)
            slv = NLOPT_LD_TNEWTON_PRECOND_RESTART!: Preconditioned truncated Newton with restarts
        case (14)
            slv = NLOPT_LD_TNEWTON_PRECOND!: Preconditioned truncated Newton
        case (15)
            slv = NLOPT_LD_TNEWTON_RESTART!: Truncated Newton with restarts
        case (16)
            slv = NLOPT_LD_TNEWTON!: Truncated Newton
        case (17)
            slv = NLOPT_LD_VAR2!: Shifted limited-memory variable-metric
        case (18)
            slv = NLOPT_LD_VAR1!: Shifted limited-memory variable-metric
        case (19)
            slv = NLOPT_LD_SLSQP!: Sequential Least-Squares Quadratic Programming (SLSQP)
        case (20)
            slv = NLOPT_GD_STOGO!: Gradient-based global optimization
        case (21)
            slv = NLOPT_GD_STOGO_RAND!: Randomized gradient-based global optimization
        case (22)
            slv = NLOPT_GN_DIRECT
        case (23)
            slv = NLOPT_GN_DIRECT_L
        case (24)
            slv = NLOPT_GN_DIRECT_L_RAND_NOSCAL
        case (25)
            slv = NLOPT_GN_DIRECT_NOSCAL
        case (26)
            slv = NLOPT_GN_DIRECT_L_NOSCAL
        case (27)
            slv = NLOPT_GN_ORIG_DIRECT
        case (28) 
            slv = NLOPT_GN_ORIG_DIRECT_L
        case (29)
            slv =  NLOPT_GN_AGS !An implementation of the algorithm AGS to solve constrained nonlinear programming problems with Lipschitzian functions.     
        case default
            slv = -100
        end select

        opt = 0

        tol = tolly
        tolc = tolly

        call nlosr(1984)

        call nlo_create(opt, slv, size_x)
        call nlo_get_algorithm(ires, opt); 
        ! call nlo_get_lower_bounds(ires, opt, lb)

        call nlo_set_population(ires, opt, 10)
       
        ! print*, "ineq  ==> ", ires
        call nlo_set_lower_bounds(ires, opt, lb)
        call nlo_get_lower_bounds(ires, opt, lb)
        call nlo_set_upper_bounds(ires, opt, ub)
        call nlo_get_upper_bounds(ires, opt, ub)
        call nlo_set_stopval(ires, opt, tol)
        call nlo_set_min_objective(ires, opt, func_, fdata)
        call nlo_set_xtol_abs(ires, opt, tol)
        call nlo_set_xtol_rel(ires, opt, tol)

        call nlo_set_maxeval(ires, opt, -1.0e0)

        call nlo_set_maxtime(ires, opt, -1.0e0)

        call nlo_set_ftol_abs(ires, opt, tol)

        call nlo_set_initial_step(ires, opt,1.0e-2*xin)
 if (present(constr)) then
 if (constr==1)call nlo_add_inequality_constraint(ires, opt, size_x, myconstraint, fdata, tolc)
 endif 
        call nlo_optimize(ires, opt, x, minf)

        call nlo_get_numevals(nevals, opt)

        call nlo_destroy(opt); 
        xin = x
        if (present(code_) .neqv. .true.) then

            print *, "Using algorithm  ==> ", ires
            print *, "pop  ==> ", ires

            print *, "lb  ==> ", lb

            print *, "ub  ==> ", ub

            print *, "stop  ==> ", ires

            write (*, '(A20,5f20.15)') ' \033[95m flag \033[0m', ires

            print *, "xtol abs  ==> ", ires

            print *, "xtol rel  ==> ", ires
            print *, "maxeval  ==> ", ires
            print *, "maxtime  ==> ", ires
            print *, "ftol  ==> ", ires
            print *, "step  ==> ", ires
            print *, "optimization  ==> ", ires
            write (*, *) '*************************************'
            write (*, '(A20,E10.3)') ' \033[95m Residual \033[0m', minf
            write (*, *) ' \033[95m Optimal x \033[0m', x
            write (*, '(A20,i)') ' \033[95m No evaluations \033[0m \033[0m', nevals
            if (ires < 0) then
                print *, ' \033[95m FAILED TO FIND OPTIMUM \033[0m'
            end if

            write (*, *) '*************************************'
        end if
    end subroutine optimize
end module
