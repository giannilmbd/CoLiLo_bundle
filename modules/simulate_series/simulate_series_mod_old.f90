module simulate_series_mod 
    use my_kinds_mod, only: dp => wp
    implicit none
    ! integer, parameter :: dp = selected_real_kind(15, 307)
    ! integer, parameter :: N = 100 ! Length of time series
    ! real(dp), dimension(:, :), allocatable :: P
    ! real(dp), dimension(:), allocatable :: X
    ! real(dp), dimension(:), allocatable :: series
    ! integer :: indx, T, i

    ! ! Example probability matrix and states vector
    ! allocate(P(2, 2))
    ! allocate(X(2))
    ! P = reshape([0.6_dp, 0.4_dp, &
    !              0.2_dp, 0.8_dp], [2, 2])
    ! X = [1.0_dp, 2.0_dp]

    ! T = 100 ! Number of time steps

    ! series = simulate_series(P, X, T)

    ! ! Print the generated time series
    ! do i = 1, N
    !     print*, series(i)
    ! end do

contains

    subroutine  simulate_series(P, X, T,series,indeces) 
        implicit none
        real(dp), dimension(:, :), intent(in) :: P
        real(dp), dimension(:), intent(in) :: X
        integer, intent(in) :: T
        real(dp), dimension(T),intent(out):: series
        integer, dimension(T),intent(out) :: indeces
        integer :: indx, i

       

        ! Choose a random initial state
        call random_number(series(1))
        indx = modulo(int(series(1) * size(X)), size(X)) + 1
        indeces(1)=indx
        series(1) = X(indx)

        do i = 2, T
            ! Get transition probabilities for the current state
            indx = modulo(indx - 1, size(X)) + 1
            ! Choose next state based on transition probabilities
            call random_number(series(i))
            indx = modulo(int(series(i) * size(X)), size(X)) + 1
            series(i) = X(indx)
            indeces(i) = indx
        end do

    end subroutine simulate_series

end module simulate_series_mod  
