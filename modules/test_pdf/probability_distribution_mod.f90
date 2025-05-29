module probability_distribution_mod
  use, intrinsic :: iso_fortran_env, only: wp => real64  
  implicit none
  real*8 :: a2
  real*8,parameter :: a1=2.0d0
  integer, parameter:: size_x=3
  integer :: i
  real*8,dimension(10) :: rngv=[(7.0*i,i=1,10)]/100.0
  real*8 :: aparams(size(rngv,1),7) ! 5 a_s, the residual and Pearson pdf at .1
 ! NOTE THAT THIS DATA IS ALSO GIVEN IN THE FUNCTION


  private i
  public calculate_residuals
  public grad
  public cfun
  interface calculate_residuals

    module procedure calculate_residuals1D
    module procedure calculate_residuals4D

end interface

interface grad
      module procedure grad1d
      module procedure grad3d
end interface


    contains
  
    ! Define your continuous PDF here
    real*8 function Pearson4( p, x)
    implicit none    
        
    real*8, intent(in) ::  p(3), x
    real*8 :: num, den, sqrt_term, atan_term, gamma_term, beta_term,a3,a4,a5,sqr,atn,gm1,gm2,bta
    logical :: flag
      
         a3=p(1)
         a4=p(2)
         a5=p(3)
      a2=a4

    flag=check_conditions( p)
    
    Pearson4=-1.
 if(flag) then
   sqr=Sqrt(-(a4**2.0/a3**2.0) + (4.0*a5)/a3)
   atn=atan((2.0*x + a4/a3)/sqr)
   gm1=Gamma((a1/(2.*a3) + 0.5*(-2.0*a2*a3 + a1*a4)/(a3**2.0*sqr)))
   gm2=Gamma(a1/(2.*a3))
   bta=Beta(-0.5 + a1/(2.*a3),0.5)
   num=(2.0**(1 - a1/a3)*(-(a4**2.0/a3**2.0) + (4.0*a5)/a3)**(-0.5 + a1/(2.*a3))*&
   exp(((-2.0*a2*a3 + a1*a4)*atn)/&
      (a3**2.0*sqr))*&
    Abs(gm1/&
       gm2)**2.0)
   den=((x**2.0 + (x*a4)/a3 + a5/a3)**(a1/(2.*a3))*bta)    
  Pearson4=num/den
 else
   print*,'\033[91m INEQUALITY CONDITIONS VIOLATED \033[0m '
 endif 
    
  end function Pearson4
  
    subroutine discretize_pdf( p,xmin, xmax, step, x_values, probabilities, n)
      real*8, intent(in),dimension(3) ::  p
      real*8 :: a3, a4, a5
      real*8, intent(in) :: xmin, xmax, step
      real*8, dimension(:), allocatable, intent(out) :: x_values, probabilities
      integer, intent(out) :: n
      real*8 :: x, pdf_sum
      integer :: i
      a3=p(1)
      a4=p(2)
      a5=p(3)

      a2=a4
      n = ceiling((xmax - xmin) / step) + 1
      allocate(x_values(n), probabilities(n))
  
      x_values = [(xmin + (i - 1) * step, i = 1, n)]
  
      probabilities(1) = Pearson4( p,xmin) * step / 2.
      probabilities(n) = Pearson4( p,xmax) * step / 2.
  
      do i = 2, n - 1
        x = x_values(i)
        probabilities(i) = Pearson4( p,x) * step
      end do
  
      ! Normalize the probabilities to ensure they sum up to 1
      pdf_sum = sum(probabilities)
      probabilities = probabilities / pdf_sum
    end subroutine discretize_pdf

    function check_conditions(p) result(flag)
        ! boundary conditions for the Pearson 4 distribution
      real*8, intent(in),dimension(3) ::  p
      real*8 :: a3, a4, a5
        logical :: flag
        real*8 :: sqtrm
        a3=p(1)
        a4=p(2)
        a5=p(3)
        a2=a4
        flag = .true.
    
        if (a3**2.0 <= 0.) then
          flag = .false.
          return
        end if
    
        if (a4**2.0 - 4. * a3 * a5 >= 0.) then
          flag = .false.
          return
        end if
    
        if (a1 / a3 <= 1.) then
          flag = .false.
          return
        end if
        sqtrm=((a2**2.0 * a3 - a1 * a2 * a4 + a1 * a4**2.0 - a3 * a4**2.0 + a1**2.0 * a5 - 4. * a1 * a3 * a5 + 4. * a3**2.0 * a5) / (a1 - 3. * a3))
        if(sqtrm<0) then
          flag = .false.
          return
        end if 
    
      end function check_conditions

subroutine cfun(m,result,n,x,grado,need_gradient,f_data)

  integer :: m,need_gradient,n
  real*8 :: result(m), x(n)
  real*8 :: grado(n,m)
  real*8 :: f_data(n),sqtrm
  real*8 ::  a3, a4, a5
  if (need_gradient.ne.0) then
    print*,'\033[36m USING PROVIDED GRADIENT \033[0m'
    grado(1,:)=[2.0*a3,4.0*a5,-a1/a3**2]
    grado(2,:)=[0.0,-2.0*a4,0.0]
    grado(3,:)=[0.0,4.0*a3,0.0]
    ! print*,grado
  endif 
  a3 =x(1)
  a4 =x(2)
  a5 =x(3)
  a2=a4
result(1)=a3**2.0 
result(2)=-(a4**2.0 - 4. * a3 * a5 )
result(3)=a1 / a3 - 1.
! sqtrm=((a2**2.0 * a3 - a1 * a2 * a4 + a1 * a4**2.0 - a3 * a4**2.0 + a1**2.0 * a5 - 4. * a1 * a3 * a5 + 4. * a3**2.0 * a5) / (a1 - 3. * a3))
! result(4)=sqtrm
print*,x
end subroutine cfun
      
        subroutine calculate_residuals1D(xin, fdata, res)
          real*8, intent(in) :: xin(size_x),fdata(size_x)
          real*8, intent(out) :: res
          
          real*8, dimension(size_x) :: res_v
          real*8 ::  a3, a4, a5 ,me,va,sk,ku
          logical :: flag
         
          ! me=fdata(1) 
          va=fdata(1)
          sk=fdata(2)
          ku=fdata(3)
          
          ! a1= xin(1)
          ! a2 =xin(2)
          a3 =xin(1)
          a4 =xin(2)
          a5 =xin(3)
          a2=a4

           flag=check_conditions( xin)
          
          ! if(flag.eqv..false.) then
          !   res_v=10000
            
          !   ! print*,'boundary not met'
          ! else 
          ! res_v(1) = me - (-a2 + a4) / (a1 - 2. * a3)
          res_v(1) = -va + (a2**2.0 * a3 - a1 * a2 * a4 + a1**2.0 * a5 + a1 * (a4**2.0 - 4. * a3 * a5) + a3 * (-a4**2.0 + 4. * a3 * a5)) / ((a1 - 3. * a3) * (a1 - 2. * a3)**2.0)
          res_v(2) = -sk + (2. * (-2. * a2 * a3 + a1 * a4) * abs(a1 - 2. * a3)) / ((a1 - 4. * a3) * (a1 - 2. * a3) * sqrt((a2**2.0 * a3 - a1 * a2 * a4 + a1 * a4**2.0 - a3 * a4**2.0 + a1**2.0 * a5 - 4. * a1 * a3 * a5 + 4. * a3**2.0 * a5) / (a1 - 3. * a3)))
          res_v(3) = -ku + (3. * (a1 - 3. * a3) * (a1**3 * a5 + 4. * a3**2.0 * (a2**2.0 + a4**2.0 - 4. * a3 * a5) - a1**2.0 * (a2 * a4 - 3. * a4**2.0 + 8. * a3 * a5) + a1 * a3 * (a2**2.0 - 4. * a2 * a4 - 5. * a4**2.0 + 20. * a3 * a5))) / ((a1 - 5. * a3) * (a1 - 4. * a3) * (a2**2.0 * a3 - a1 * a2 * a4 + a1**2.0 * a5 + a1 * (a4**2.0 - 4. * a3 * a5) + a3 * (-a4**2.0 + 4. * a3 * a5)))
            
        ! end if
       
          res = sqrt(dot_product(res_v, res_v))

             if(flag.eqv..false.) then
            ! res=abs(exp(4*res)-1.0)*1000.0
               res=9999.99
            ! print*,'boundary not met'
             endif
            
            ! 
! print*,xin
          
        end subroutine calculate_residuals1D


        subroutine calculate_residuals4D(xin, fdata, res_v)
          real*8, intent(in) :: xin(size_x),fdata(size_x)
          real*8, dimension(size_x) :: res_v
          real*8 ::  a3, a4, a5,me,va,sk,ku
          logical :: flag

          
          va=fdata(1)
          sk=fdata(2)
          ku=fdata(3)
          
          ! a1= xin(1)
          ! a2 =xin(2)
          a3 =xin(1)
          a4 =xin(2)
          a5 =xin(3)
          a2=a4
          flag=check_conditions( xin)
          
          if(flag.eqv..false.) then
            res_v=10000
            
            ! print*,'boundary not met'
          else 
          ! res_v(1) = me - (-a2 + a4) / (a1 - 2. * a3)
          res_v(1) = -va + (a2**2.0 * a3 - a1 * a2 * a4 + a1**2.0 * a5 + a1 * (a4**2.0 - 4. * a3 * a5) + a3 * (-a4**2.0 + 4. * a3 * a5)) / ((a1 - 3. * a3) * (a1 - 2. * a3)**2.0)
          res_v(2) = -sk + (2. * (-2. * a2 * a3 + a1 * a4) * abs(a1 - 2. * a3)) / ((a1 - 4. * a3) * (a1 - 2. * a3) * sqrt((a2**2.0 * a3 - a1 * a2 * a4 + a1 * a4**2.0 - a3 * a4**2.0 + a1**2.0 * a5 - 4. * a1 * a3 * a5 + 4. * a3**2.0 * a5) / (a1 - 3. * a3)))
          res_v(3) = -ku + (3. * (a1 - 3. * a3) * (a1**3 * a5 + 4. * a3**2.0 * (a2**2.0 + a4**2.0 - 4. * a3 * a5) - a1**2.0 * (a2 * a4 - 3. * a4**2.0 + 8. * a3 * a5) + a1 * a3 * (a2**2.0 - 4. * a2 * a4 - 5. * a4**2.0 + 20. * a3 * a5))) / ((a1 - 5. * a3) * (a1 - 4. * a3) * (a2**2.0 * a3 - a1 * a2 * a4 + a1**2.0 * a5 + a1 * (a4**2.0 - 4. * a3 * a5) + a3 * (-a4**2.0 + 4. * a3 * a5)))
           
        end if
       
        end subroutine calculate_residuals4D        
      
        subroutine min_resid(result, n, x, grado, need_gradient, f_data)
          integer :: need_gradient,n
          real*8 :: result, x(n)
          real*8 :: grado(n)
          real*8 :: f_data(n)
          if (need_gradient.ne.0) then
            print*,'\033[36m USING PROVIDED GRADIENT \033[0m'
            call grad(x,grado)
            print*,grado
          endif 
          call calculate_residuals(x,f_data,result)
       
          end  subroutine min_resid

          function f4fzero(x,f_data)
            ! call fzero(x_both, foc2, check,parin=parin0)
            real*8 , intent(in):: x(:)
            real*8 , intent(in):: f_data(:)
            real*8 :: f4fzero(size(x,1))
            
            
            call calculate_residuals(x,f_data,f4fzero)
            
          end  function f4fzero
 
  pure real*8 function beta(a, b)
  real*8, intent(in) :: a, b

  ! Compute the Beta function using the gamma function from the standard library
  beta = (gamma(a)* gamma(b))/gamma(a + b)

end function beta

subroutine grad3d(xin,grad)
  implicit none
  real*8 :: grad(size_x,size_x)
  real*8 :: vgrad(size_x*size_x)
  real*8 :: a3,a4,a5
  real*8, intent(in),dimension(size_x) :: xin
  
  a3 =xin(1)
  a4 =xin(2)
  a5 =xin(3)
  a2 = a4;

  vgrad=[(3d0*a5)/(a1 - 3d0*a3)**2,0.d0,1/(a1 - 3d0*a3),((5d0*a1 - 12d0*a3)*a4*Sqrt(a5/(a1 - 3d0*a3)))/((a1 - 4d0*a3)**2d0*a5),&
  2/((a1 - 4d0*a3)*Sqrt(a5/(a1 - 3d0*a3))),-(a4/((a1 - 4d0*a3)*a5*Sqrt(a5/(a1 - 3d0*a3)))),&
  (6*(60d0*a3**2d0*a4**2 + a1**3d0*a5 + a1**2*(6d0*a4**2 - 8d0*a3*a5) + 8d0*a1*a3*(-5d0*a4**2 + 2d0*a3*a5)))/&
   ((a1 - 5d0*a3)**2*(a1 - 4d0*a3)**2d0*a5),(12*(a1 - 3d0*a3)*a4)/((a1 - 5d0*a3)*(a1 - 4d0*a3)*a5),&
  (-6*(a1 - 3d0*a3)*a4**2)/((a1 - 5d0*a3)*(a1 - 4d0*a3)*a5**2)]
  grad=Transpose(reshape(source=vgrad,shape=shape(grad)))
end subroutine grad3d

subroutine grad1d(xin,grad)
  implicit none
  real*8,intent(out) :: grad(size_x)
  
  real*8 :: a3,a4,a5
  real*8, intent(in),dimension(size_x) :: xin
  
  a3 =xin(1)
  a4 =xin(2)
  a5 =xin(3)
  a2 = a4;

  grad=[(a1**12*a5*(18.*a5 - 5.*a4*Sqrt(a5/(a1 - 3*a3))) + &
     a1**11*(154.*a4**2*a5 - 576.*a3*a5**2 + 192.*a3*a4*a5*Sqrt(a5/(a1 - 3*a3))) + &
     a1**2*a3**8*(2.0030976e7*a4**4 - 5.792256e7*a3*a4**2*a5 + 6.096384e6*a3**2*a5**2 + 1.90368e6*a5**4 - &
        2.140344e7*a3**2*a4*a5*Sqrt(a5/(a1 - 3*a3))) + &
     a3**10*(2.79936e6*a4**4 - 6.89472e6*a3*a4**2*a5 + 383999.9999999999*a5**4 - &
        2.592e6*a3**2*a4*a5*Sqrt(a5/(a1 - 3*a3))) + &
     a1**10*(216.*a4**4 - 4664.*a3*a4**2*a5 + 8315.999999999998*a3**2*a5**2 + 3.*a5**4 - &
        3336.9999999999995*a3**2*a4*a5*Sqrt(a5/(a1 - 3*a3))) + &
     a1*a3**9*(-1.119744e7*a4**4 + 2.9746656000000004e7*a3*a4**2*a5 - 1.492992e6*a3**2*a5**2 - 1.2864e6*a5**4 + &
        1.10592e7*a3**2*a4*a5*Sqrt(a5/(a1 - 3*a3))) + &
     a1**4*(-5.1560504e7*a3**7*a4**2*a5 + a3**6*(1.4487984e7*a4**4 + 907728.*a5**4) + &
        a3**8*a5*(1.2312864e7*a5 - 1.9241032e7*a4*Sqrt(a5/(a1 - 3*a3)))) + &
     a1**6*(-1.0389672e7*a3**5*a4**2*a5 + a3**4*(2.1846959999999995e6*a4**4 + 86052.*a5**4) + &
        a3**6*a5*(4.499856e6*a5 - 4.113971e6*a4*Sqrt(a5/(a1 - 3*a3)))) + &
     a1**8*(-517672.*a3**3*a4**2*a5 + a3**2*(68544.*a4**4 + 1629.*a5**4) + &
        a3**4*a5*(406674.*a5 - 240711.*a4*Sqrt(a5/(a1 - 3*a3)))) + &
     a1**9*(63690.*a3**2*a4**2*a5 + a3*(-5760.*a4**4 - 105.*a5**4) + &
        a3**3*a5*(-71496.*a5 + 34712.*a4*Sqrt(a5/(a1 - 3*a3)))) + &
     a1**7*(2.782862e6*a3**4*a4**2*a5 + a3**3*(-479520.*a4**4 - 14739.*a5**4) + &
        a3**5*a5*(-1.606824e6*a5 + 1.172464e6*a4*Sqrt(a5/(a1 - 3*a3)))) + &
     a1**5*(2.7491326e7*a3**6*a4**2*a5 + a3**5*(-6.7752e6*a4**4 - 338472.*a5**4) + &
        a3**7*a5*(-8.931456e6*a5 + 1.047996e7*a4*Sqrt(a5/(a1 - 3*a3)))) + &
     a1**3*(6.7180656e7*a3**8*a4**2*a5 + a3**7*(-2.109888e7*a4**4 - 1.638192e6*a5**4) + &
        a3**9*a5*(-1.1228544e7*a5 + 2.4837216e7*a4*Sqrt(a5/(a1 - 3*a3)))))/&
   ((a1 - 5.*a3)**3*(a1 - 4.*a3)**3*(a1 - 3.*a3)**3*(a1 - 2.*a3)**4*a5**2*&
     Sqrt(a5**2/(a1 - 3*a3)**2 + (9*(a1 - 3*a3)**2*(2*a4**2 + (a1 - 4*a3)*a5)**2)/&
        ((a1 - 5*a3)**2*(a1 - 4*a3)**2*a5**2) + (1 - (2*a4)/((a1 - 4*a3)*Sqrt(a5/(a1 - 3*a3))))**2)),&
  (2591.9999999999995*a3**4*a4**3 - 6383.999999999999*a3**5*a4*a5 - 2.*a1**6*a5*Sqrt(a5/(a1 - 3*a3)) - &
     2400.*a3**6*a5*Sqrt(a5/(a1 - 3*a3)) + a1**5*(40.*a4*a5 + 42.*a3*a5*Sqrt(a5/(a1 - 3*a3))) + &
     a1**4*(72.*a4**3 - 572.*a3*a4*a5 - 357.99999999999994*a3**2*a5*Sqrt(a5/(a1 - 3*a3))) + &
     a1**3*(-720.*a3*a4**3 + 3216.*a3**2*a4*a5 + 1581.9999999999998*a3**3*a5*Sqrt(a5/(a1 - 3*a3))) + &
     a1**2*(2663.9999999999995*a3**2*a4**3 - 8876.*a3**3*a4*a5 - 3816.*a3**4*a5*Sqrt(a5/(a1 - 3*a3))) + &
     a1*(-4320.*a3**3*a4**3 + 12015.999999999998*a3**4*a4*a5 + 4760.*a3**5*a5*Sqrt(a5/(a1 - 3*a3))))/&
   ((a1 - 5.*a3)**2*(a1 - 4.*a3)**2*(a1 - 2.*a3)**2*a5**2*&
     Sqrt(a5**2/(a1 - 3*a3)**2 + (9*(a1 - 3*a3)**2*(2*a4**2 + (a1 - 4*a3)*a5)**2)/&
        ((a1 - 5*a3)**2*(a1 - 4*a3)**2*a5**2) + (1 - (2*a4)/((a1 - 4*a3)*Sqrt(a5/(a1 - 3*a3))))**2)),&
  (1.*a1**10*a4*a5*Sqrt(a5/(a1 - 3*a3)) + a1**9*(-20.*a4**2*a5 - 31.*a3*a4*a5*Sqrt(a5/(a1 - 3*a3))) + &
     a1**8*(-36.*a4**4 + 486.*a3*a4**2*a5 + 1.*a5**4 + 425.99999999999994*a3**2*a4*a5*Sqrt(a5/(a1 - 3*a3))) + &
     a3**8*(-46656.*a4**4 + 114912.*a3*a4**2*a5 + 6400.*a5**4 + 43200.*a3**2*a4*a5*Sqrt(a5/(a1 - 3*a3))) + &
     a1**7*(720.*a3*a4**4 - 5207.999999999999*a3**2*a4**2*a5 - 26.*a3*a5**4 - &
        3417.9999999999995*a3**3*a4*a5*Sqrt(a5/(a1 - 3*a3))) + &
     a1*(-407808.*a3**8*a4**2*a5 - 157680.*a3**9*a4*a5*Sqrt(a5/(a1 - 3*a3)) + &
        a3**7*(155520.*a4**4 - 18560.*a5**4)) + &
     a1**3*(-578384.*a3**6*a4**2*a5 - 243015.99999999997*a3**7*a4*a5*Sqrt(a5/(a1 - 3*a3)) + &
        a3**5*(185760.*a4**4 - 16000.*a5**4)) + &
     a1**5*(-127763.99999999999*a3**4*a4**2*a5 - 62222.99999999999*a3**5*a4*a5*Sqrt(a5/(a1 - 3*a3)) + &
        a3**3*(30960.*a4**4 - 1791.9999999999998*a5**4)) + &
     a1**6*(32300.*a3**3*a4**2*a5 + 17736.999999999996*a3**4*a4*a5*Sqrt(a5/(a1 - 3*a3)) + &
        a3**2*(-6263.999999999999*a4**4 + 288.99999999999994*a5**4)) + &
     a1**4*(334254.*a3**5*a4**2*a5 + 149500.*a3**6*a4*a5*Sqrt(a5/(a1 - 3*a3)) + &
        a3**4*(-95076.*a4**4 + 6776.*a5**4)) + &
     a1**2*(638351.9999999999*a3**7*a4**2*a5 + 255887.99999999997*a3**8*a4*a5*Sqrt(a5/(a1 - 3*a3)) + &
        a3**6*(-225504.*a4**4 + 23055.999999999996*a5**4)))/&
   ((a1 - 5.*a3)**2*(a1 - 4.*a3)**2*(a1 - 3.*a3)**2*(a1 - 2.*a3)**4*a5**3*&
     Sqrt(a5**2/(a1 - 3*a3)**2 + (9*(a1 - 3*a3)**2*(2*a4**2 + (a1 - 4*a3)*a5)**2)/&
        ((a1 - 5*a3)**2*(a1 - 4*a3)**2*a5**2) + (1 - (2*a4)/((a1 - 4*a3)*Sqrt(a5/(a1 - 3*a3))))**2))]

end subroutine grad1d
  end module probability_distribution_mod
  