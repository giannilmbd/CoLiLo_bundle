module MarkovMoments_mod
    use my_kinds_mod, only: wp
    use globals
    implicit none
    
 

    
contains

subroutine stationary_distribution(P, pi)
    real(kind=wp), intent(in) :: P(NS,NS)
    real(kind=wp), intent(out) :: pi(NS)
    real(kind=wp) :: A(NS+1,NS), e(NS+1),identity_N(NS,NS)
    integer :: info,cnt

    identity_N=0.0_wp
    do cnt=1,NS
        identity_N(cnt,cnt)=1.0_wp
    enddo 
    ! Construct the matrix A for solving P'*pi= pi and sum(pi) = 1
    A(1:NS, 1:NS) = (identity_N- transpose(P))
    A(NS+1, 1:NS) = 1.0_wp ! ROW NS+1 is of ones
    ! print*,A(1:NS,1:NS)
    ! print*,A(NS+1,:)
    e = 0.0_wp
    e(NS+1) = 1.0_wp

    ! Solve the system of equations to find the stationary distribution
    call solve_system(A, e, pi)
    ! if((sum(pi,1)>1.0e-6_wp).or.(sum(pi,1)<(0.999999))) then
    !     print*, ' ERROR: ERGODIC PROBABILITY SUMS .ne. 1'
    !     print*, 'Sum of probabilities',sum(pi,1)
    !     error stop
    ! endif 
    
end subroutine stationary_distribution

    
subroutine solve_system(A, b, x)
    use my_kinds_mod, only: wp, sp, dp
    implicit none
    real(kind=wp), intent(in) :: A(NS+1,NS), b(NS+1)
    real(kind=wp), intent(out) :: x(NS)
    integer :: ipiv(NS), info, Lwork
    real(kind=wp), allocatable, dimension(:) :: work

    ! Solve the system of linear equations Ax = b
    Lwork = -1
    allocate(work(NS))

    if (wp == sp) then
        call sgels('N', NS+1, NS, 1, A, NS+1, b, NS+1, work, Lwork, info)
    else if (wp == dp) then
        call dgels('N', NS+1, NS, 1, A, NS+1, b, NS+1, work, Lwork, info)
    else
        print *, "Unsupported precision for wp"
        stop
    end if

    Lwork = work(1)
    deallocate(work)

    allocate(work(Lwork))

    if (wp == sp) then
        call sgels('N', NS+1, NS, 1, A, NS+1, b, NS+1, work, Lwork, info)
    else if (wp == dp) then
        call dgels('N', NS+1, NS, 1, A, NS+1, b, NS+1, work, Lwork, info)
    else
        print *, "Unsupported precision for wp"
        stop
    end if

    deallocate(work)
    x = b(1:NS)
    if (info /= 0) then
        print *, "Error solving system of equations"
        stop
    end if
end subroutine solve_system
    
    subroutine compute_moment(P, X, k,moment)
        real(kind=wp), intent(in) :: P(NS,NS), X(NS)
        real(kind=wp),  intent(in) :: k
        real(kind=wp),intent(out) :: moment
        real(kind=wp) :: pi(NS)
        integer :: i
        
        ! Compute the kth order uncentered moment
        call stationary_distribution(P,pi)
        ! moment = 0.0
        ! do i = 1, NS
        !     moment = moment + pi(i) * X(i)**k
        ! end do

        moment=dot_product(pi,X**k)
        
    end subroutine compute_moment

end module MarkovMoments_mod
