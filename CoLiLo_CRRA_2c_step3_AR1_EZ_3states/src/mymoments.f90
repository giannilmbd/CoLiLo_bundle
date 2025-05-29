module mymoments
    ! This modules generates sample mean, stdev, variance, skewness and kurtosis
    use my_kinds_mod, only: wp
    implicit none
    contains 
   pure function mean(x)
        implicit none
        real(kind=wp) :: mean
        real(kind=wp), intent(in)::x(:)
        real(kind=wp) :: T
    
        T=size(x)
        mean=sum(x)/T
    end function
    pure function stdev(x)
        implicit none
        real(kind=wp) :: stdev
        real(kind=wp), intent(in)::x(:)

        stdev=sqrt(variance(x))
    end function
    pure function cov(x,y)
        real(kind=wp) :: cov
        real(kind=wp), intent(in) :: x(:),y(:)
        integer :: T

        T=size(x)
        cov=dot_product((x-mean(x)),(y-mean(y)))/T
       


    end function
    pure function variance(x)
        implicit none
        real(kind=wp) :: variance
        real(kind=wp), intent(in)::x(:)
        integer :: T

        T=size(x)
        variance=sum((x-mean(x))**2)/T
    end function
    pure function skewness(x)
        implicit none
        real(kind=wp) :: skewness,mean_x,stdev_x
        real(kind=wp), intent(in)::x(:)
        integer :: T
    
        T=size(x)
        skewness=(sum((x-mean(x))**3)/T)/(stdev(x)**3)
    end function
    pure function kurtosis(x)
        implicit none
        real(kind=wp) :: kurtosis
        real(kind=wp), intent(in)::x(:)
        integer :: T
    
        T=size(x)
        kurtosis=(sum((x-mean(x))**4)/T)/(stdev(x)**4)
    end function
end module mymoments