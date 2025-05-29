module probability_distribution_mod
  use, intrinsic :: iso_fortran_env, only: wp => real64  
  implicit none
    contains
  
    ! Define your continuous PDF here
    real(wp) function Pearson4(a1, a2, a3, a4, a5, x)
        
        
    real(wp), intent(in) :: a1, a2, a3, a4, a5, x
    real(wp) :: num, den, sqrt_term, atan_term, gamma_term, beta_term
    logical :: flag
    flag=check_conditions(a1, a2, a3, a4, a5)
    Pearson4=-1._wp
 if(flag) then
  Pearson4=(2.0**(1 - a1/a3)*(-(a4**2.0/a3**2) + (4.0*a5)/a3)**(-0.5 + a1/(2.*a3))*&
      exp(((-2.0*a2*a3 + a1*a4)*atan((2.0*x + a4/a3)/Sqrt(-(a4**2.0/a3**2.0) + (4.0*a5)/a3)))/&
         (a3**2.0*Sqrt(-(a4**2.0/a3**2.0) + (4.0*a5)/a3)))*&
       Abs(Gamma((a1/(2.*a3) + 0.5*(-2.0*a2*a3 + a1*a4)/(a3**2*Sqrt(-(a4**2.0/a3**2.0) + (4.0*a5)/a3))))/&
          Gamma(a1/(2.*a3)))**2)/((x**2 + (x*a4)/a3 + a5/a3)**(a1/(2.*a3))*Beta(-0.5 + a1/(2.*a3),0.5_wp))
 endif 
    
  end function Pearson4
  
    subroutine discretize_pdf(a1, a2, a3, a4, a5,xmin, xmax, step, x_values, probabilities, n)
      real(wp), intent(in) :: a1, a2, a3, a4, a5
      real(wp), intent(in) :: xmin, xmax, step
      real(wp), dimension(:), allocatable, intent(out) :: x_values, probabilities
      integer, intent(out) :: n
  
      real(wp) :: x, pdf_sum
      integer :: i
  
      n = ceiling((xmax - xmin) / step) + 1
      allocate(x_values(n), probabilities(n))
  
      x_values = [(xmin + (i - 1) * step, i = 1, n)]
  
      probabilities(1) = Pearson4(a1, a2, a3, a4, a5,xmin) * step / 2._wp
      probabilities(n) = Pearson4(a1, a2, a3, a4, a5,xmax) * step / 2._wp
  
      do i = 2, n - 1
        x = x_values(i)
        probabilities(i) = Pearson4(a1, a2, a3, a4, a5,x) * step
      end do
  
      ! Normalize the probabilities to ensure they sum up to 1
      pdf_sum = sum(probabilities)
      probabilities = probabilities / pdf_sum
    end subroutine discretize_pdf

    function check_conditions(a1, a2, a3, a4, a5) result(flag)
        ! boundary conditions for the Pearson 4 distribution
        real(wp), intent(in) :: a1, a2, a3, a4, a5
        logical :: flag
    
        flag = .true.
    
        if (a3**2 <= 0._wp) then
          flag = .false.
          return
        end if
    
        if (a4**2 - 4._wp * a3 * a5 >= 0._wp) then
          flag = .false.
          return
        end if
    
        if (a1 / a3 <= 1._wp) then
          flag = .false.
          return
        end if
    
      end function check_conditions


      
        subroutine calculate_residuals(xin, fdata, res)
          real(wp), intent(in) :: xin(5),fdata(4)
          real(wp), intent(out) :: res
          real(wp), dimension(4) :: res_v
          real(wp) :: a1, a2, a3, a4, a5,me,va,sk,ku
          logical :: flag
          me=fdata(1) 
          va=fdata(2)
          sk=fdata(3)
          ku=fdata(4)
          a1= xin(1)
          a2 =xin(2)
          a3 =xin(3)
          a4 =xin(4)
          a5 =xin(5)

          flag=check_conditions(a1, a2, a3, a4, a5)
          
      
          res_v(1) = me - (-a2 + a4) / (a1 - 2._wp * a3)
          res_v(2) = -va + (a2**2 * a3 - a1 * a2 * a4 + a1**2 * a5 + a1 * (a4**2 - 4._wp * a3 * a5) + a3 * (-a4**2 + 4._wp * a3 * a5)) / ((a1 - 3._wp * a3) * (a1 - 2._wp * a3)**2)
          res_v(3) = -sk + (2._wp * (-2._wp * a2 * a3 + a1 * a4) * abs(a1 - 2._wp * a3)) / ((a1 - 4._wp * a3) * (a1 - 2._wp * a3) * sqrt((a2**2 * a3 - a1 * a2 * a4 + a1 * a4**2 - a3 * a4**2 + a1**2 * a5 - 4._wp * a1 * a3 * a5 + 4._wp * a3**2 * a5) / (a1 - 3._wp * a3)))
          res_v(4) = -ku + (3._wp * (a1 - 3._wp * a3) * (a1**3 * a5 + 4._wp * a3**2 * (a2**2 + a4**2 - 4._wp * a3 * a5) - a1**2 * (a2 * a4 - 3._wp * a4**2 + 8._wp * a3 * a5) + a1 * a3 * (a2**2 - 4._wp * a2 * a4 - 5._wp * a4**2 + 20._wp * a3 * a5))) / ((a1 - 5._wp * a3) * (a1 - 4._wp * a3) * (a2**2 * a3 - a1 * a2 * a4 + a1**2 * a5 + a1 * (a4**2 - 4._wp * a3 * a5) + a3 * (-a4**2 + 4._wp * a3 * a5)))
          res = sqrt(dot_product(res_v, res_v))
          if(flag.eqv..false.) then
            res=res+10000
            print*,'boundary not met'
          endif 
          
        end subroutine calculate_residuals
      
        subroutine min_resid(result, n, x, grad, need_gradient, f_data)
          integer :: need_gradient,n
          real(wp) :: result, x(n)
          real(wp) :: grad(n)
          real(wp) :: f_data(4)
          if (need_gradient.ne.0) then
            print*,"No Gradient computed"
          endif 
          call calculate_residuals(x,f_data,result)
       
          end  subroutine min_resid

 
  pure real(wp) function beta(a, b)
  real(wp), intent(in) :: a, b

  ! Compute the Beta function using the gamma function from the standard library
  beta = (gamma(a)* gamma(b))/gamma(a + b)

end function beta

  end module probability_distribution_mod
  