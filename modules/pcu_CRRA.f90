module pcu_CRRA_module
    use globals
    use my_kinds_mod
    implicit none
contains
    pure function pcu_CRRA(W, EC, EL) result(pcu)
! compute permanent consumption units as
! W=welfare under CM - Welfare under autarky
! Find pcu such that W=U_AU(pcu*C_AU,L_AU)-U_AU(C_AU,L_AU)
        real(kind=wp), intent(in) :: W, EC, EL
        real(kind=wp) :: pcu

! W= (pc1*EC1)**(1.0_wp-1.0_wp/gamma)/(1.0_wp-1.0_wp/gamma)-chhi*EL1**(1.0_wp+varphi)/(1.0_wp+varphi)
! (pc1*EC1)**(1.0_wp-1.0_wp/gamma)/(1.0_wp-1.0_wp/gamma)=W+chhi*EL1**(1.0_wp+varphi)/(1.0_wp+varphi)
! pcu=(((W+chhi*EL**(1.0_wp+varphi)/(1.0_wp+varphi))*(1.0_wp-1.0_wp/gamma))**(1/(1.0_wp-1.0_wp/gamma)))/EC

        pcu = (((W + EC**(1.0_wp - 1.0_wp/gamma)/(1.0_wp - 1.0_wp/gamma))*(1.0_wp - 1.0_wp/gamma))**(1/(1.0_wp - 1.0_wp/gamma)))/EC

    end function
end module
