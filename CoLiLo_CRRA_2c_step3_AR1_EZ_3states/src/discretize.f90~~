
module discretize
private :: myinit_random_seed
  ! use ifport ! this was used for random number, now  standard Fortran and ce toolbox (for seed)
!#####################################################
!#from the codes produced by the authors of “Discretizing a Process with Non-zero Skewness and High Kurtosis”, 2015
!# 		Simone Civale, Luis Díez-Catalán, Fatih Fazilet
!#
!#####################################################            
    CONTAINS

    SUBROUTINE TransitionOneStep(inputpar,states,transition)
      !# state transition matrix
      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(9)                          ::  inputpar
      REAL*8, DIMENSION(:), INTENT(IN)                          ::  states
      REAL*8, DIMENSION(size(states),size(states)), INTENT(OUT) ::  transition
      !
      REAL*8                                                    ::  m1,s1,m2,s2,m3,s3,p1,p2,rho
      REAL*8, DIMENSION(SIZE(states)-1)                         ::  nodes
      INTEGER                                                     ::  i,j,n
      !
      m1=inputpar(1)          ! mean1
      s1=inputpar(2)**0.5d0  ! sd1
      m2=inputpar(3)          ! mean2
      s2=inputpar(4)**0.5d0  ! sd2
      m3=inputpar(5)          ! mean3    
      s3=inputpar(6)**0.5d0  ! sd3
      p1=inputpar(7)          ! p1
      p2=inputpar(8)          ! p2
      rho=inputpar(9)         ! rho
      transition=0
      n=SIZE(states)
      DO i=1,n-1
          nodes(i)=(states(i+1)+states(i))/2.0d0
      END DO
      DO i=1,n
          DO j=1,n
              IF (j .EQ. 1) THEN
                  transition(i,j)=alnorm((nodes(1)-rho*states(i)-m1)/s1,.false.)*p1 +      &
                                & alnorm((nodes(1)-rho*states(i)-m2)/s2,.false.)*p2 +      &
                                & alnorm((nodes(1)-rho*states(i)-m3)/s3,.false.)*(1.0d0-p1-p2)
                  ELSE IF (j .EQ. n) THEN
                      transition(i,j)=alnorm((nodes(n-1)-rho*states(i)-m1)/s1,.true.)*p1 + &
                                    & alnorm((nodes(n-1)-rho*states(i)-m2)/s2,.true.)*p2 + &
                                    & alnorm((nodes(n-1)-rho*states(i)-m3)/s3,.true.)*(1.0d0-p1-p2)    
                  ELSE
                      transition(i,j)=alnorm((nodes(j)   - rho*states(i)-m1)/s1,.false.)*p1 - &
                                    & alnorm((nodes(j-1) - rho*states(i)-m1)/s1,.false.)*p1 + &
                                    & alnorm((nodes(j)   - rho*states(i)-m2)/s2,.false.)*p2 - &
                                    & alnorm((nodes(j-1) - rho*states(i)-m2)/s2,.false.)*p2 + &
                                    & alnorm((nodes(j)   - rho*states(i)-m3)/s3,.false.)*(1.0d0-p1-p2) - &
                                    & alnorm((nodes(j-1) - rho*states(i)-m3)/s3,.false.)*(1.0d0-p1-p2)
              END IF                                              
          END DO
      END DO
!         print*,transition(1,:)           
!         print*,transition(2,:)           
!         print*,transition(3,:)
!         print*,states
!         print*,''      
  END SUBROUTINE TransitionOneStep

function myrand(x)
  implicit none
  integer,optional :: x
  real*8 :: myrand

  if(present(x)) then
    call myinit_random_seed(.true.) 
  end if
  call random_number(myrand)
end function myrand

  SUBROUTINE simulate(CumTransition,SimulatedProcess)
      IMPLICIT NONE
      REAL*8, DIMENSION(:,:), INTENT(IN)  ::  CumTransition     ! the cumulated transition matrix
      REAL*8, DIMENSION(:), INTENT(OUT)   ::  SimulatedProcess  ! the realizations of the continuous process
      ! this subroutine returns a vector of simulated data from a markov process characterized by states and transition probabilities
      INTEGER                                     ::  n                 ! number of states used for the discrete approximation
      INTEGER                                     ::  t1                ! number of simulations
      INTEGER                                     ::  i               ! counters
      REAL*8, DIMENSION(SIZE(CumTransition,1))  ::  difference        ! housekeeping
      REAL*8                                    ::  prob
      !
      
      prob=myrand(3456)
      n=SIZE(CumTransition,1)
      t1=SIZE(SimulatedProcess)
      ! we set the first simulated value
      SimulatedProcess(1)=FLOOR(n/2.0d0)
      ! we make the values fluctuate according to the transition probabilities
      DO i=2,t1
          prob=myrand()
          difference=CumTransition(INT(SimulatedProcess(i-1)),:)-prob
          SimulatedProcess(i)=MINLOC(difference,1,difference .GT. 0.0d0)
      END DO
  END SUBROUTINE simulate    
function alnorm ( x, upper )
    !*****************************************************************************80
    !
    !! ALNORM computes the cumulative density of the standard normal distribution.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    13 January 2008
    !
    !  Author:
    !
    !    Original FORTRAN77 version by David Hill.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    David Hill,
    !    Algorithm AS 66:
    !    The Normal Integral,
    !    Applied Statistics,
    !    Volume 22, Number 3, 1973, pages 424-427.
    !
    !  Parameters:
    !
    !    Input, REAL*8 X, is one endpoint of the semi-infinite interval
    !    over which the integration takes place.
    !
    !    Input, logical UPPER, determines whether the upper or lower
    !    interval is to be integrated:
    !    .TRUE.  => integrate from X to + Infinity;
    !    .FALSE. => integrate from - Infinity to X.
    !
    !    Output, REAL*8 ALNORM, the integral of the standard normal
    !    distribution over the desired interval.
    !
      implicit none

      REAL*8, parameter :: a1 = 5.75885480458D+00
      REAL*8, parameter :: a2 = 2.62433121679D+00
      REAL*8, parameter :: a3 = 5.92885724438D+00
      REAL*8 alnorm
      REAL*8, parameter :: b1 = -29.8213557807D+00
      REAL*8, parameter :: b2 = 48.6959930692D+00
      REAL*8, parameter :: c1 = -0.000000038052D+00
      REAL*8, parameter :: c2 = 0.000398064794D+00
      REAL*8, parameter :: c3 = -0.151679116635D+00
      REAL*8, parameter :: c4 = 4.8385912808D+00
      REAL*8, parameter :: c5 = 0.742380924027D+00
      REAL*8, parameter :: c6 = 3.99019417011D+00
      REAL*8, parameter :: con = 1.28D+00
      REAL*8, parameter :: d1 = 1.00000615302D+00
      REAL*8, parameter :: d2 = 1.98615381364D+00
      REAL*8, parameter :: d3 = 5.29330324926D+00
      REAL*8, parameter :: d4 = -15.1508972451D+00
      REAL*8, parameter :: d5 = 30.789933034D+00
      REAL*8, parameter :: ltone = 7.0D+00
      REAL*8, parameter :: p = 0.398942280444D+00
      REAL*8, parameter :: q = 0.39990348504D+00
      REAL*8, parameter :: r = 0.398942280385D+00
      logical up
      logical upper
      REAL*8, parameter :: utzero = 18.66D+00
      REAL*8 x
      REAL*8 y
      REAL*8 z

      up = upper
      z = x

      if ( z < 0.0D+00 ) then
        up = .not. up
        z = - z
      end if

      if ( ltone < z .and. ( ( .not. up ) .or. utzero < z ) ) then

        if ( up ) then
          alnorm = 0.0D+00
        else
          alnorm = 1.0D+00
        end if

        return

      end if

      y = 0.5D+00 * z * z

      if ( z <= con ) then

        alnorm = 0.5D+00 - z * ( p - q * y &
          / ( y + a1 + b1 &
          / ( y + a2 + b2 &
          / ( y + a3 ))))

      else

        alnorm = r * exp ( - y ) &
          / ( z + c1 + d1 &
          / ( z + c2 + d2 &
          / ( z + c3 + d3 &
          / ( z + c4 + d4 &
          / ( z + c5 + d5 &
          / ( z + c6 ))))))

      end if

      if ( .not. up ) then
        alnorm = 1.0D+00 - alnorm
      end if

      return
    end


    
    FUNCTION cumulative(x)
        REAL*8, DIMENSION(:,:), INTENT(IN)      :: x
        REAL*8, DIMENSION(size(x,1),size(x,2))  :: cumulative
        ! this function returns the cumulative of a transition matrix x
        INTEGER                                   :: row,col
        INTEGER                                   :: counter1,counter2
        row=size(x,1)
        col=size(x,2)
        DO counter1=1,row
            DO counter2=1,col
                IF (counter2 .EQ. 1) THEN
                    cumulative(counter1,counter2)=x(counter1,counter2)
                ELSE
                    cumulative(counter1,counter2)=cumulative(counter1,counter2-1)+x(counter1,counter2)
                END IF
            END DO
        END DO
    END FUNCTION            
    FUNCTION linspace(a,b,num)
      ! generates a linearly spaced support for the distribution
        IMPLICIT NONE
        REAL*8, INTENT(IN)                :: a,b
        INTEGER, INTENT(IN)                 :: num
        REAL*8, DIMENSION(num)             :: linspace
        real*8                    ::  rng,stp,tmp(num)
        ! linspace generates a linearly spaced grid in the interval [a,b]
        ! with the number of points being determined by num. Clearly a < b is required and
        ! num > 0. The output is x, the linearly spaced grid
        INTEGER  :: i
        ! rng=b-a
        ! stp=rng/(num-1)
        ! tmp=a
        ! do i=2,(num-1)
        !   tmp(i)=tmp(i-1)+stp
        ! enddo 
        ! tmp(num)=b
        ! linspace=tmp

        linspace = (/ (a+(real(i,8)-1.0d0)*b/(real(num,8)-1.0d0), i = 1, num) /)
    END FUNCTION linspace


    FUNCTION randnorm()
      IMPLICIT NONE
      REAL*8    :: randnorm
      ! returns a normally distributed number with zero mean and unit variance
      REAL*8    :: rsq,v1,v2
      !
      DO
          v1=myrand()
          v2=myrand()
          v1=2.0d0*v1-1.0d0
          v2=2.0d0*v2-1.0d0
          rsq=v1**2+v2**2
          IF (rsq > 0.0 .and. rsq < 1.0) EXIT
      END DO
      rsq=sqrt(-2.0d0*log(rsq)/rsq)
      randnorm=v1*rsq
  END FUNCTION randnorm

  subroutine myinit_random_seed(fixed)
     
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, dt(8)
    integer, parameter :: int64 = selected_int_kind(16)
    integer(kind=int64) :: t
    logical, optional :: fixed
 
    call random_seed(size = n)
    allocate(seed(n))
 
    call system_clock(t)
    if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
            + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
            + dt(3) * 24_int64 * 60 * 60 * 1000 &
            + dt(5) * 60 * 60 * 1000 &
            + dt(6) * 60 * 1000 + dt(7) * 1000 &
            + dt(8)
    endif
    
    if(present(fixed))then
        if(fixed)t = 0
    endif
    do i = 1, n
        seed(i) = lcg(t)
    enddo
    call random_seed(put=seed)
 
contains
 
    ! This simple PRNG might not be good enough for real work, but is
    ! sufficient for seeding a better PRNG.
    function lcg(s)
 
        implicit none
        integer :: lcg
        integer(int64) :: s
 
        if (s == 0) then
            s = 104729
        else
            s = mod(s, 4294967296_int64)
        endif
        s = mod(s * 279470273_int64, 4294967291_int64)
        lcg = int(mod(s, int(huge(0), int64)), kind(0))
 
    end function lcg
 
end subroutine myinit_random_seed

    end module discretize    