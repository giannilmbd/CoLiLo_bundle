module pcu_CRRA_module
    use globals
    use my_kinds_mod
    implicit none
contains
    pure function pcu_CRRA(W, EC) result(pcu)
! compute permanent consumption units as
! W=welfare under CM - Welfare under autarky
! Find pcu such that W=U_AU(pcu*C_AU,L_AU)-U_AU(C_AU,L_AU)
        real(kind=wp), intent(in) :: W, EC
        real(kind=wp) :: pcu

! W= (pc1*EC1)**(1.0_wp-1.0_wp/gamma)/(1.0_wp-1.0_wp/gamma)-(EC1)**(1.0_wp-1.0_wp/gamma)/(1.0_wp-1.0_wp/gamma)
! pc1=((((W*(1.0_wp-1.0_wp/gamma) EC**( 1.0_wp/gamma-1.0_wp)+1)**(1/(1.0_wp-1.0_wp/gamma))

        pcu = ((W *(1.0_wp - 1.0_wp/gamma)*EC**( 1.0_wp/gamma-1.0_wp))+1.0_wp)**(1.0_wp/(1.0_wp - 1.0_wp/gamma))

    end function
end module
