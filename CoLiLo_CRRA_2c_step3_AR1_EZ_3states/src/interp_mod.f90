module interp_mod
    use bspline_sub_module       ! provides db4ink, db4val
    use bspline_kinds_module, only: ip
    use globals
    implicit none
    private
    public :: bracket, bilinear_eval, quadlinear_eval
    public :: bspline4_init , bspline4_eval
  
  contains
  
    !---------------------------------------------------------------------
    integer function bracket(grid, n, x)
      real(kind=wp), intent(in) :: grid(:)
      integer,   intent(in)    :: n
      real(kind=wp), intent(in) :: x
      integer :: lo, hi, mid
      if (x <= grid(1)) then
        bracket = 1
      elseif (x >= grid(n)) then
        bracket = n - 1
      else
        lo = 1; hi = n
        do while (hi - lo > 1)
          mid = (lo + hi) / 2
          if (x < grid(mid)) then
            hi = mid
          else
            lo = mid
          end if
        end do
        bracket = lo
      end if
    end function bracket
  
    !---------------------------------------------------------------------
    real(kind=wp) function bilinear_eval(x, y, xg, ng, yg, mg, f)
      real(kind=wp), intent(in) :: x, y
      real(kind=wp), intent(in) :: xg(ng), yg(mg)
      integer,   intent(in) :: ng, mg
      real(kind=wp), intent(in) :: f(ng, mg)
      integer :: i, j
      real(kind=wp) :: tx, ty
      i = bracket(xg, ng, x)
      j = bracket(yg, mg, y)
      tx = (x - xg(i)) / (xg(i+1) - xg(i))
      ty = (y - yg(j)) / (yg(j+1) - yg(j))
      bilinear_eval = (1.0_wp - tx)*(1.0_wp - ty)*f(i  , j  ) + &
                      tx        *(1.0_wp - ty)*f(i+1, j  ) + &
                      (1.0_wp - tx)*       ty      *f(i  , j+1) + &
                      tx        *       ty      *f(i+1, j+1)
    end function bilinear_eval
  
    !---------------------------------------------------------------------
    real(kind=wp) function quadlinear_eval(x1, x2, x3, x4, &
                                          g1, n1, g2, n2, g3, n3, g4, n4, f4) result(val)
      real(kind=wp), intent(in) :: x1, x2, x3, x4
      real(kind=wp), intent(in) :: g1(n1), g2(n2), g3(n3), g4(n4)
      integer, intent(in) :: n1, n2, n3, n4
      real(kind=wp), intent(in) :: f4(n1, n2, n3, n4)
      integer :: b1, b2, b3, b4, c1, c2, c3, c4
      real(kind=wp) :: w1lo, w1hi, w2lo, w2hi, w3lo, w3hi, w4lo, w4hi
  
      b1 = bracket(g1, n1, x1)
      b2 = bracket(g2, n2, x2)
      b3 = bracket(g3, n3, x3)
      b4 = bracket(g4, n4, x4)
  
      w1hi = (x1 - g1(b1)) / (g1(b1+1) - g1(b1)); w1lo = 1.0_wp - w1hi
      w2hi = (x2 - g2(b2)) / (g2(b2+1) - g2(b2)); w2lo = 1.0_wp - w2hi
      w3hi = (x3 - g3(b3)) / (g3(b3+1) - g3(b3)); w3lo = 1.0_wp - w3hi
      w4hi = (x4 - g4(b4)) / (g4(b4+1) - g4(b4)); w4lo = 1.0_wp - w4hi
  
      val = 0.0_wp
      do c1 = 0, 1
        do c2 = 0, 1
          do c3 = 0, 1
            do c4 = 0, 1
              val = val + merge(w1lo, w1hi, c1==1) * &
                        merge(w2lo, w2hi, c2==1) * &
                        merge(w3lo, w3hi, c3==1) * &
                        merge(w4lo, w4hi, c4==1) * &
                        f4(b1+c1, b2+c2, b3+c3, b4+c4)
            end do
          end do
        end do
      end do
    end function quadlinear_eval
  
    !---------------------------------------------------------------------
    subroutine bspline4_init(grid1, grid2, grid3, grid4, values, &
        tx,ty,tz,tq, bcoef4, iflag)
        use bspline_kinds_module, only: wp,ip
        use bspline_module, only: db4ink
        implicit none
      
        real(wp), intent(in) :: grid1(0:NK1), grid2(0:NK2), grid3(0:NK3), grid4(0:NK4)
        real(wp), intent(inout) :: values(0:NK1,0:NK2,0:NK3,0:NK4)

        real(wp), intent(out) :: bcoef4(size(grid1), size(grid2), size(grid3), size(grid4))
        integer(ip), intent(out)  :: iflag
        integer:: deg_int2
        integer(ip):: n1,n2,n3,n4,iknot
        integer ::  max_knot_len
        real(wp) :: tx(NK1+1+deg_int), ty(NK2+1+deg_int), tz(NK3+1+deg_int), tq(NK4+1+deg_int)
      
        n1 = INT(size(grid1),ip)
        n2 = INT(size(grid2),ip)
        n3 = INT(size(grid3),ip)
        n4 = INT(size(grid4),ip)
        iknot = 0 ! automatically select the knots
        deg_int2 = INT(deg_int,ip)
        call db4ink(grid1, n1, grid2, n2, grid3, n3, grid4, n4, values, &
                    deg_int2, deg_int2, deg_int2, deg_int2, iknot, tx, ty, tz, tq, bcoef4, iflag)
                    ! print*,"IFLAG",iflag
        ! if(iflag.ne.INT(0,ip)) then
        !     print*,"deg_int2=",deg_int2
        !     print*,"IFLAG",iflag
        !   write(*,*) ANSI_RED//'Error in bspline4_init: '//get_status_message(iflag)//ANSI_RESET
        ! !   stop
        ! end if
      
      end subroutine bspline4_init
      
    !---------------------------------------------------------------------
      real(wp) function bspline4_eval(x_state, tx, ty, tz, tq, bcoef4,  extrap) result(val)
      use bspline_kinds_module, only: wp, ip
      use bspline_module, only: db4val
      implicit none
    
      real(wp), intent(in) :: x_state(4)
      real(wp) :: tx(NK1+1+deg_int), ty(NK2+1+deg_int), tz(NK3+1+deg_int), tq(NK4+1+deg_int)
      
      real(wp), intent(inout) :: bcoef4(:,:,:,:)
      integer(ip) :: k
      logical, intent(in)  :: extrap
    
      integer(ip) :: idx, idy, idz, idq, valflag
      integer(ip) :: inbvx, inbvy, inbvz, inbvq
      integer(ip) :: iloy, iloz, iloq
      integer(ip) :: n1, n2, n3, n4
      real(wp), allocatable :: w3(:), w2(:,:), w1(:,:,:), w0(:)
      real(wp) :: ftmp
      k = deg_int
      n1 = size(bcoef4, 1)
      n2 = size(bcoef4, 2)
      n3 = size(bcoef4, 3)
      n4 = size(bcoef4, 4)
    
      allocate(w3(k), w2(k,k), w1(k,k,k), w0(3*k))
    
      idx = 0; idy = 0; idz = 0; idq = 0
      valflag = 0
      inbvx = 1; inbvy = 1; inbvz = 1; inbvq = 1
      iloy = 1; iloz = 1; iloq = 1
    
      call db4val(x_state(1), x_state(2), x_state(3), x_state(4), &
                  idx, idy, idz, idq, &
                  tx, ty, tz, tq, n1, n2, n3, n4, &
                  k, k, k, k, bcoef4, ftmp, valflag, &
                  inbvx, inbvy, inbvz, inbvq, &
                  iloy, iloz, iloq, &
                  w3, w2, w1, w0, extrap)
      val = ftmp
    
      deallocate(w3, w2, w1, w0)
    end function bspline4_eval
  
  end module interp_mod
  