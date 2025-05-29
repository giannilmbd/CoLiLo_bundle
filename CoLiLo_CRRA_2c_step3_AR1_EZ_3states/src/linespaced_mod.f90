module linespaced_mod
    use my_kinds_mod, only: wp
    implicit none

    contains

    function linespaced_values(Nin, LV, UV) result(values)
        implicit none
        integer, intent(in) :: Nin
        real(kind=wp), intent(in) :: LV, UV
        real(kind=wp), dimension(Nin) :: values
        real(kind=wp) :: delta
        integer :: i

        delta = (UV - LV) / (Nin  - 1.0_wp)
        values(1)=LV
        values(Nin)=UV
        do i = 2, Nin-1 
            values(i) = LV + (real(i,wp) - 1.0_wp) * delta
        end do

    end function linespaced_values

end module linespaced_mod