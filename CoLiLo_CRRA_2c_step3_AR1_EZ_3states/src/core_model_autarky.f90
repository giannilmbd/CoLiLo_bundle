subroutine core_model_autarky(pf, C, Cs, L, Ls, ph, Q, Y1, res, pred, D, Ds)
    ! implicit none
    use my_kinds_mod, only: wp
    use globals, dump=>core_model_autarky

    implicit none

    real(kind=wp), intent(in)  ::   D, Ds
    real(kind=wp)              :: pf, a, b
    integer             :: cntr
    real(kind=wp), intent(out) :: C, Cs, L, Ls, ph, Q, pred, res, Y1

    res = 1.0e10_wp

    if (trade_elast < 100.0e0_wp) then
        ph = ((1.0e0_wp - (1.0e0_wp - nu)*pf**(1.0e0_wp - trade_elast))**(1.0e0_wp/(1.0e0_wp - trade_elast)))*(nu)**(1.0e0_wp/(trade_elast-1.0e0_wp ))

        Q = ((1.0e0_wp-nu_b)*ph**(1.0e0_wp - trade_elast) + nu_b*pf**(1.0e0_wp - trade_elast))**(1.0e0_wp/(1.0e0_wp - trade_elast)); 
        
        a = (1.0e0_wp - alphha)/(varphi + alphha + (1.0e0_wp - alphha)/gamma)
        
        b = (varphi + 1)/(varphi + alphha + (1.0e0_wp - alphha)/gamma)

        if (alphha == 1.0e0_wp) then

            C = D*ph
            Cs = Ds*pf/Q
            L = 1.0e0_wp

            Ls = 1.0e0_wp
        else
            C = ((1.0e0_wp - alphha)/chhi)**a*D**b*ph**b
            Cs = ((1.0e0_wp - alphha)/chhi)**a*Ds**b*pf**b*Q**(-b)

            L = C**((-1.0e0_wp/(gamma*(varphi + alphha))))*(((1.0e0_wp - alphha)/chhi)*ph*(D))**(1.0e0_wp/(varphi + alphha))

            Ls = Cs**((-1.0e0_wp/(gamma*(varphi + alphha))))*(((1.0e0_wp - alphha)/chhi)*pf/Q*(Ds))**(1.0e0_wp/(varphi + alphha))
        end if

        Y1 = (D)*L**(1.0e0_wp - alphha)
        !

    pred = (((1.0e0_wp - nu)*n/(1.0e0_wp - n)*C + nu_b*Q**(trade_elast)*Cs)/((Ds)*Ls**(1.0e0_wp - alphha)))**(1.0e0_wp/trade_elast) !pf**(1.0e0_wp-trade_elast)*
    else

        pred = 1.0e0_wp
        ph = pred
        Q = 1.0e0_wp
        a = (1.0e0_wp - alphha)/(varphi + alphha + (1.0e0_wp - alphha)/gamma)
        b = (varphi + 1)/(varphi + alphha + (1.0e0_wp - alphha)/gamma)

        if (alphha == 1.0e0_wp) then

            C = D*ph
            Cs = Ds*pf/Q
            L = 1.0e0_wp

            Ls = 1.0e0_wp
        else
            C = ((1.0e0_wp - alphha)/chhi)**a*D**b*ph**b
            Cs = ((1.0e0_wp - alphha)/chhi)**a*Ds**b*pf**b*(Q)**(-b)

            L = C**((-1.0e0_wp/(gamma*(varphi + alphha))))*(((1.0e0_wp - alphha)/chhi)*ph*(D))**(1.0e0_wp/(varphi + alphha))

            Ls = Cs**((-1.0e0_wp/(gamma*(varphi + alphha))))*(((1.0e0_wp - alphha)/chhi)*pred/Q*(Ds))**(1.0e0_wp/(varphi + alphha))
        end if

        Y1 = (D)*L**(1.0e0_wp - alphha)
    end if
    res = pred/pf - 1.0e0_wp

end subroutine core_model_autarky


