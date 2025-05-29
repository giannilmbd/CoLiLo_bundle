subroutine welfare(welf,kappa)

    use my_toolbox
    use globals
    implicit none

    real(kind=wp), dimension(4), intent(out) :: welf
    real(kind=wp), intent(in) :: kappa


    ! you might want to average over possible initial conditions or use the same initial condition as for the Arrow sec

    !welf(1) = ((1.0d0 - 1.0d0/gamma)*spline_eval([1.0d0, 1.0d0, 1.0d0], &
                          ! coeff_v1(:, :, :, pos_zero1, pos_zero2), LB_all, UB_all))**(1/((1.0d0 - 1.0d0/gamma)))
    !welf(2) = ((1.0d0 - 1.0d0/gamma)*spline_eval([1.0d0, 1.0d0, 1.0d0], &
                          ! coeff_v2(:, :, :, pos_zero1, pos_zero2), LB_all, UB_all))**(1/((1.0d0 - 1.0d0/gamma)))
													
													
	! Note that welfare is already scaled correctly: V1 AND V2 are the pdv of utility.				
    welf(1) = -spline_eval([A1(pos_zeroA1),A2(pos_zeroA2),OmegaS(pos_zeroOm)], coeff_v1(:,:,:, pos_zero1,pos_zero2), LB_all, UB_all) 																
    welf(2) = -spline_eval([A1(pos_zeroA1),A2(pos_zeroA2),OmegaS(pos_zeroOm)], coeff_v2(:,:,:, pos_zero1,pos_zero2), LB_all, UB_all) 	

! autarky
    ! call value_autarky()

  !  welf(3) = ((1.0d0 - 1.0d0/gamma)*spline_eval([1.0d0,1.0d0], &
   !                       coeff_v1_aut(:, :, pos_zero1, pos_zero2), LB_all(2:3), UB_all(2:3)))**(1/(1.0d0 - 1.0d0/gamma))
    !welf(4) = ((1.0d0 - 1.0d0/gamma)*spline_eval([1.0d0,1.0d0], &
     !                     coeff_v2_aut(:, :, pos_zero1, pos_zero2), LB_all(2:3), UB_all(2:3)))**(1/((1.0d0 - 1.0d0/gamma)))
  welf(3)= -v1_aut(pos_zeroA1,pos_zeroA2,pos_zero1,pos_zero2)												
  welf(4)= -v2_aut(pos_zeroA1,pos_zeroA2,pos_zero1,pos_zero2)
end subroutine
