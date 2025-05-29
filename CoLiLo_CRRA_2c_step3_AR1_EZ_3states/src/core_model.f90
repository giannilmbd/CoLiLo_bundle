subroutine core_model(pf, kappa, OmegaS_in, C1, C2, L1, L2, ph, Q, Y1, res, pred, D1, D2)
    ! implicit none
    use globals, dump => core_model
    implicit none
    ! integer, intent(in) :: rows,cols
    ! real(kind=wp), intent(in), dimension(rows,cols) :: x
    ! integer, intent(in) :: timezero
    real(kind=wp), intent(in):: kappa, pf,   D1, D2,OmegaS_in
    real(kind=wp)            :: mu, mus, YW,temp1,temp2 ,exp1,exp2 ,term1,term2
    integer           :: cntr
    real(kind=wp), intent(out) :: C1, C2, L1, L2, ph, Q, pred, res, Y1

    ! real(kind=wp), dimension(0:ngrid_-1) :: predX, predp
    res = 1.0d10
! print*,'nu',nu
! print*,'pf',pf
! read(*,*) Y1
! cntr=0.0e0_wp
!  iter_alloc: do while ((abs(res)>sig).and.(cntr<itermax))
!  cntr=cntr+1.0e0_wp
    ! COBB-DOUGLAS CASE
    if (trade_elast == 1.0e0_wp) then
        ph = pf**((nu - 1.0e0_wp)/nu)
        Q = ph**((1.0e0_wp-nu_b))*pf**(nu_b)
        mu=1.0e0_wp/(n+(1.0e0_wp-n)*Q**(1-gamma)*kappa*OmegaS_in**gamma)
        mus =  mu*kappa*Q**(-gamma)*OmegaS_in
        if (alphha == 1.0e0_wp) then
            YW = n*ph*D1 + (1.0e0_wp - n)*pf*D2
        else
            YW = (n*ph*D1*(mu**(-1.0e0_wp/gamma)*(1.0e0_wp - alphha)/chhi*ph*D1)**((1.0e0_wp - alphha)/(varphi + alphha)) + &
        (1.0e0_wp-n)*pf*D2*(mus**(-1.0e0_wp/gamma)*(1.0e0_wp -alphha )/chhi*pf/Q*D2)**((1.0e0_wp-alphha)/(varphi+alphha)))**((varphi+alphha)/((varphi+alphha)+(1-alphha)/gamma))
        end if

        C1 = mu*YW
        C2 = mus*YW

        if (alphha == 1.0e0_wp) then
            ! L1=log(exp(C1)**((-rho)/(varphi+alphha))*((1-alphha)/chhi*exp(ph)*exp(D1))**(1/(varphi+alphha)));
            L1 = 1.0e0_wp

            ! L2=log(exp(C2)**((-rho)/(varphi+alphha))*((1-alphha)/chhi*exp(pf)/exp(Q)*exp(D2))**(1/(varphi+alphha)));
            L2 = 1.0e0_wp
        else
            ! L1=log(exp(C1)**((-rho)/(varphi+alphha))*((1-alphha)/chhi*exp(ph)*exp(D1))**(1/(varphi+alphha)));
            L1 = C1**((-1.0e0_wp/(gamma*(varphi + alphha))))*(((1.0e0_wp - alphha)/chhi)*ph*(D1))**(1.0e0_wp/(varphi + alphha))

            ! L2=log(exp(C2)**((-rho)/(varphi+alphha))*((1.0e0_wp-alphha)/chhi*exp(pf)/exp(Q)*exp(D2))**(1.0e0_wp/(varphi+alphha)));
            L2 = C2**((-1.0e0_wp/(gamma*(varphi + alphha))))*(((1.0e0_wp - alphha)/chhi)*pf/Q*(D2))**(1.0e0_wp/(varphi + alphha))
        end if

        Y1 = (D1)*L1**(1.0e0_wp - alphha)
        ! Welf1 = exp(C1)**(1 - rho)/(1 - rho) - chhi*exp(L1)**(1 + varphi)/(1 + varphi)
        ! Welf2 = exp(C2)**(1 - rho)/(1 - rho) - chhi*exp(L2)**(1 + varphi)/(1 + varphi)
        !######################################################  COMPUTE RESIDUALS
        !     predX=(XX+1)*((exp(C1)**(1-rho)-exp(C1)**(-rho)*exp(ph)*exp(D1)*exp(L1)**(1-alphha)+betta*EXX)**2+1)/(XX**2+1)-1
        ! res(0)=DNRM2(ngrid_,(XX**2+1)/(predX**2+1)-1,1) !

        pred = (((1.0e0_wp - nu)*n/(1.0e0_wp - n)*C1 + (nu_b)*Q*C2)/((D2)*L2**(1.0e0_wp - alphha))) !pf**(1.0e0_wp-trade_elast)*
! CES CASE
    elseif ((trade_elast.ne.1.0e0_wp).and.(trade_elast < 100.0e0_wp)) then

        ph = ((1.0e0_wp - (1.0e0_wp - nu)*pf**(1.0e0_wp - trade_elast))/nu)**(1.0e0_wp/(1.0e0_wp - trade_elast))

        Q = ((1.0e0_wp-nu_b)*ph**(1.0e0_wp - trade_elast) + (nu_b)*pf**(1.0e0_wp - trade_elast))**(1.0e0_wp/(1.0e0_wp - trade_elast)); 
        
        ! mu=1.0e0_wp/(n*(nu*ph**(1-trade_elast)+(1-nu)*pf**(1-trade_elast))+&
        ! (1.0e0_wp-n)*((1.0e0_wp-nu_b)*ph**(1-trade_elast)+nu_b*pf**(1-trade_elast))*Q**(trade_elast-gamma)*kappa)
        mu=1.0e0_wp/(n+(1.0e0_wp-n)*Q**(1-gamma)*kappa*OmegaS_in**gamma)
        
        mus = mu*kappa*Q**(-gamma)*OmegaS_in
        
        if (alphha == 1.0e0_wp) then
            YW = n*ph*D1 + (1.0e0_wp - n)*pf*D2
        else

            temp1 = mu**(-1.0e0_wp / gamma) * (1.0e0_wp - alphha) / chhi * ph * D1
            temp2 = mus**(-1.0e0_wp / gamma) * (1.0e0_wp - alphha) / chhi * pf / Q * D2

            term1 = (n * ph * D1 * temp1**((1.0e0_wp - alphha) / (varphi + alphha)))
            term2 = ((1.0e0_wp - n) * pf * D2 * temp2**((1.0e0_wp - alphha) / (varphi + alphha)))

            YW = (term1 + term2)**((varphi + alphha) / ((varphi + alphha) + (1.0e0_wp - alphha) / gamma))

        !     YW = (n*ph*D1*(mu**(-1.0e0_wp/gamma)*(1.0e0_wp - alphha)/chhi*ph*D1)**((1.0e0_wp - alphha)/(varphi + alphha)) + &
        ! (1.0e0_wp-n)*pf*D2*(mus**(-1.0e0_wp/gamma)*(1.0e0_wp - alphha)/chhi*pf/Q*D2)**((1.0e0_wp-alphha)/(varphi+alphha)))**((varphi+alphha)/((varphi+alphha)+(1-alphha)/gamma))
        end if

        C1  = mu*YW

        C2 = mus*YW

        if (alphha == 1.0e0_wp) then
            ! L1=log(exp(C1)**((-rho)/(varphi+alphha))*((1-alphha)/chhi*exp(ph)*exp(D1))**(1/(varphi+alphha)));
            L1 = 1.0e0_wp

            ! L2=log(exp(C2)**((-rho)/(varphi+alphha))*((1-alphha)/chhi*exp(pf)/exp(Q)*exp(D2))**(1/(varphi+alphha)));
            L2 = 1.0e0_wp
        else

            exp1 = -1.0e0_wp / (gamma * (varphi + alphha))
            exp2 = 1.0e0_wp / (varphi + alphha)

            L1 = C1**exp1 * (((1.0e0_wp - alphha) / chhi) * ph * D1)**exp2
            L2 = C2**exp1 * (((1.0e0_wp - alphha) / chhi) * pf / Q * D2)**exp2
            ! L1=log(exp(C1)**((-rho)/(varphi+alphha))*((1-alphha)/chhi*exp(ph)*exp(D1))**(1/(varphi+alphha)));
            ! L1 = C1**((-1.0e0_wp/(gamma*(varphi + alphha))))*(((1.0e0_wp - alphha)/chhi)*ph*(D1))**(1.0e0_wp/(varphi + alphha))

            ! ! L2=log(exp(C2)**((-rho)/(varphi+alphha))*((1.0e0_wp-alphha)/chhi*exp(pf)/exp(Q)*exp(D2))**(1.0e0_wp/(varphi+alphha)));
            ! L2 = C2**((-1.0e0_wp/(gamma*(varphi + alphha))))*(((1.0e0_wp - alphha)/chhi)*pf/Q*(D2))**(1.0e0_wp/(varphi + alphha))
        end if

        Y1 = (D1)*L1**(1.0e0_wp - alphha)
        ! Welf1 = exp(C1)**(1 - rho)/(1 - rho) - chhi*exp(L1)**(1 + varphi)/(1 + varphi)
        ! Welf2 = exp(C2)**(1 - rho)/(1 - rho) - chhi*exp(L2)**(1 + varphi)/(1 + varphi)
        !######################################################  COMPUTE RESIDUALS
        !     predX=(XX+1)*((exp(C1)**(1-rho)-exp(C1)**(-rho)*exp(ph)*exp(D1)*exp(L1)**(1-alphha)+betta*EXX)**2+1)/(XX**2+1)-1
        ! res(0)=DNRM2(ngrid_,(XX**2+1)/(predX**2+1)-1,1) !

pred = (((1.0e0_wp - nu)*n/(1.0e0_wp - n)*C1 + (nu_b)*Q**(trade_elast)*C2)/((D2)*L2**(1.0e0_wp - alphha)))**(1.0e0_wp/trade_elast) !pf**(1.0e0_wp-trade_elast)*
    ! PERFECT SUBSTITUTION    
    else

        pred = 1.0e0_wp
        ph = pred
        Q = 1.0e0_wp
        mu = 1.0e0_wp/(n + (1.0e0_wp - n)*kappa*OmegaS_in**gamma)
        mus =  mu*kappa*Q**(-gamma)*OmegaS_in
        if (alphha == 1.0e0_wp) then
            YW = n*ph*D1 + (1.0e0_wp - n)*pred*D2
        else
            YW = (n*ph*D1*(mu**(-1.0e0_wp/gamma)*(1.0e0_wp - alphha)/chhi*ph*D1)**((1.0e0_wp - alphha)/(varphi + alphha)) + &
        (1.0e0_wp-n)*pred*D2*(mus**(-1.0e0_wp/gamma)*(1.0e0_wp -alphha )/chhi*pred/Q*D2)**((1.0e0_wp-alphha)/(varphi+alphha)))**((varphi+alphha)/((varphi+alphha)+(1-alphha)/gamma))
        end if
        C1 = mu*YW
        C2 = mus*YW

        if (alphha == 1.0e0_wp) then
            ! L1=log(exp(C1)**((-rho)/(varphi+alphha))*((1-alphha)/chhi*exp(ph)*exp(D1))**(1/(varphi+alphha)));
            L1 = 1.0e0_wp

            ! L2=log(exp(C2)**((-rho)/(varphi+alphha))*((1-alphha)/chhi*exp(pf)/exp(Q)*exp(D2))**(1/(varphi+alphha)));
            L2 = 1.0e0_wp
        else
            ! L1=log(exp(C1)**((-rho)/(varphi+alphha))*((1-alphha)/chhi*exp(ph)*exp(D1))**(1/(varphi+alphha)));
            L1 = C1**((-1.0e0_wp/(gamma*(varphi + alphha))))*(((1.0e0_wp - alphha)/chhi)*ph*(D1))**(1.0e0_wp/(varphi + alphha))

            ! L2=log(exp(C2)**((-rho)/(varphi+alphha))*((1.0e0_wp-alphha)/chhi*exp(pf)/exp(Q)*exp(D2))**(1.0e0_wp/(varphi+alphha)));
            L2 = C2**((-1.0e0_wp/(gamma*(varphi + alphha))))*(((1.0e0_wp - alphha)/chhi)*pred/Q*(D2))**(1.0e0_wp/(varphi + alphha))
        end if

        Y1 = (D1)*L1**(1.0e0_wp - alphha)
    end if
    res = 2.0_wp * (pred / pf - 1.0e0_wp) / (pred / pf + 1.0e0_wp)

!         pf=smooth_pf*pred+(1.0e0_wp-smooth_pf)*pf
!  print*,'pf core',pf
!      end do iter_alloc

end subroutine core_model


