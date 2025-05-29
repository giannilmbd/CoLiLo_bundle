module my_spline_mod
   use my_kinds_mod, only: wp
   implicit none

   private
   public :: spline_interp, spline_eval, spline

   interface spline_interp
       module procedure spline_interp1, spline_interp2, spline_interp3, spline_interp4
   end interface

   interface spline_eval
       module procedure &
           spline1, spline1_grid, &
           spline2, spline2_grid, &
           spline3, spline3_grid, &
           spline4, spline4_grid, &
           spline5, spline5_grid, &
           spline6, spline6_grid, &
           spline7, spline7_grid
   end interface

   interface spline
       module procedure &
           spline1_complete, spline1_complete_m, &
           spline2_complete, spline2_complete_m, &
           spline3_complete, spline3_complete_m, &
           spline4_complete, spline4_complete_m, &
           spline5_complete, spline5_complete_m, &
           spline6_complete, spline6_complete_m, &
           spline7_complete, spline7_complete_m
   end interface

contains

!===============================================================================
! Helper: assert_eq
!===============================================================================
function assert_eq(n1, n2, string) result(res)
   integer, intent(in) :: n1, n2
   character(len=*), intent(in) :: string
   integer :: res

   if (n1 == n2) then
       res = n1
   else
       call error(string, 'assert_eq: assertion failed')
   end if
end function assert_eq

!===============================================================================
! Helper: error
!===============================================================================
subroutine error(string1, string2)
   character(len=*), intent(in) :: string1, string2
   write(*,*) 'Error in routine:', trim(string1)
   write(*,*) 'Message:', trim(string2)
   stop
end subroutine error
!===============================================================================
! 1D Spline Interpolation
!===============================================================================
subroutine spline_interp1(yi, c)
   real(wp), intent(in)  :: yi(0:)
   real(wp), intent(out) :: c(1:)
   real(wp) :: r(assert_eq(size(yi,1)+2, size(c,1), 'spline_interp1'))
   real(wp) :: d(assert_eq(size(yi,1)+2, size(c,1), 'spline_interp1'))
   integer :: n, j

   n = assert_eq(size(yi,1)+2, size(c,1), 'spline_interp1')
   n = n - 3

   ! Numerical derivatives at endpoints
   r(1)    = (2d0*yi(0)-5d0*yi(1)+4d0*yi(2)-yi(3))/6d0
   r(n+3)  = (2d0*yi(n)-5d0*yi(n-1)+4d0*yi(n-2)-yi(n-3))/6d0

   r(2:n+2) = yi(0:n)

   ! Solve tridiagonal system
   c(2) = (yi(0)-r(1))/6d0
   c(n+2) = (yi(n)-r(n+3))/6d0

   d(3) = 4d0
   r(3) = yi(1) - c(2)
   r(n+1) = yi(n-1) - c(n+2)

   do j = 4, n+1
       d(j) = 4d0 - 1d0/d(j-1)
       r(j) = r(j) - r(j-1)/d(j-1)
   end do

   c(n+1) = r(n+1)/d(n+1)
   do j = n, 3, -1
       c(j) = (r(j) - c(j+1))/d(j)
   end do

   c(1) = r(1) + 2d0*c(2) - c(3)
   c(n+3) = r(n+3) + 2d0*c(n+2) - c(n+1)
end subroutine spline_interp1

!===============================================================================
! 2D Spline Interpolation
!===============================================================================
subroutine spline_interp2(yi, c)
   real(wp), intent(in)  :: yi(:,:)
   real(wp), intent(out) :: c(:,:)
   integer :: i

   do i = 1, size(yi,1)
       call spline_interp1(yi(i,:), c(i,:))
   end do
end subroutine spline_interp2

!===============================================================================
! 3D Spline Interpolation
!===============================================================================
subroutine spline_interp3(yi, c)
   real(wp), intent(in)  :: yi(:,:,:)
   real(wp), intent(out) :: c(:,:,:)
   integer :: i, j

   do i = 1, size(yi,1)
       do j = 1, size(yi,2)
           call spline_interp1(yi(i,j,:), c(i,j,:))
       end do
   end do
end subroutine spline_interp3

!===============================================================================
! 4D Spline Interpolation
!===============================================================================
subroutine spline_interp4(yi, c)
   real(wp), intent(in)  :: yi(:,:,:,:)
   real(wp), intent(out) :: c(:,:,:,:)
   integer :: i, j, k

   do i = 1, size(yi,1)
       do j = 1, size(yi,2)
           do k = 1, size(yi,3)
               call spline_interp1(yi(i,j,k,:), c(i,j,k,:))
           end do
       end do
   end do
end subroutine spline_interp4
!===============================================================================
! 1D Spline Evaluation
!===============================================================================
function spline1(x, c, lb, ub) result(y_interp)
   real(wp), intent(in) :: x(:)
   real(wp), intent(in) :: c(1:)
   real(wp), intent(in) :: lb(:), ub(:)
   real(wp) :: y_interp

   integer :: n, ix
   real(wp) :: t

   n = size(c) - 3

   if (x(1) <= 0.0_wp) then
       y_interp = c(2) + (c(3)-c(2))*x(1)
   else if (x(1) >= real(n-1,wp)) then
       y_interp = c(n+1) + (c(n+2)-c(n+1))*(x(1)-real(n-1,wp))
   else
       ix = int(x(1))
       t = x(1) - real(ix,wp)
       y_interp = ( &
           (c(ix+1)*(1.0_wp-t)**3 + &
            c(ix+2)*(3.0_wp*t**3 - 6.0_wp*t**2 + 4.0_wp) + &
            c(ix+3)*(-3.0_wp*t**3 + 3.0_wp*t**2 + 3.0_wp*t + 1.0_wp) + &
            c(ix+4)*t**3) / 6.0_wp )
   end if
end function spline1

function spline1_grid(x, c, lb, ub) result(y_interp)
   real(wp), intent(in) :: x(:)
   real(wp), intent(in) :: c(1:)
   real(wp), intent(in) :: lb(:), ub(:)
   real(wp) :: y_interp

   y_interp = spline1(x, c, lb, ub)
end function spline1_grid

!===============================================================================
! 2D Spline Evaluation
!===============================================================================
function spline2(x, c, lb, ub) result(y_interp)
   real(wp), intent(in) :: x(:)
   real(wp), intent(in) :: c(:,:)
   real(wp), intent(in) :: lb(:), ub(:)
   real(wp) :: y_interp

   integer :: ix1
   real(wp) :: dx1
   real(wp) :: c0(1:size(c,2))

   ix1 = int(x(1))
   dx1 = x(1) - real(ix1,wp)

   c0 = (1.0_wp-dx1)*c(ix1,:) + dx1*c(ix1+1,:)
   y_interp = spline1([x(2)], c0, lb, ub)
end function spline2

function spline2_grid(x, c, lb, ub) result(y_interp)
   real(wp), intent(in) :: x(:)
   real(wp), intent(in) :: c(:,:)
   real(wp), intent(in) :: lb(:), ub(:)
   real(wp) :: y_interp

   y_interp = spline2(x, c, lb, ub)
end function spline2_grid

!===============================================================================
! 3D Spline Evaluation
!===============================================================================
function spline3(x, c, lb, ub) result(y_interp)
   real(wp), intent(in) :: x(:)
   real(wp), intent(in) :: c(:,:,:)
   real(wp), intent(in) :: lb(:), ub(:)
   real(wp) :: y_interp

   integer :: ix1, ix2
   real(wp) :: dx1, dx2
   real(wp) :: c00(1:size(c,3)), c01(1:size(c,3))
   real(wp) :: y0, y1

   ix1 = int(x(1))
   dx1 = x(1) - real(ix1,wp)

   ix2 = int(x(2))
   dx2 = x(2) - real(ix2,wp)

   c00 = (1.0_wp-dx1)*c(ix1,ix2,:) + dx1*c(ix1+1,ix2,:)
   c01 = (1.0_wp-dx1)*c(ix1,ix2+1,:) + dx1*c(ix1+1,ix2+1,:)

   y0 = spline1([x(3)], c00, lb, ub)
   y1 = spline1([x(3)], c01, lb, ub)

   y_interp = (1.0_wp-dx2)*y0 + dx2*y1
end function spline3

function spline3_grid(x, c, lb, ub) result(y_interp)
   real(wp), intent(in) :: x(:)
   real(wp), intent(in) :: c(:,:,:)
   real(wp), intent(in) :: lb(:), ub(:)
   real(wp) :: y_interp

   y_interp = spline3(x, c, lb, ub)
end function spline3_grid

!===============================================================================
! 4D Spline Evaluation
!===============================================================================
function spline4(x, c, lb, ub) result(y_interp)
   real(wp), intent(in) :: x(:)
   real(wp), intent(in) :: c(:,:,:,:)
   real(wp), intent(in) :: lb(:), ub(:)
   real(wp) :: y_interp

   integer :: ix1, ix2, ix3
   real(wp) :: dx1, dx2, dx3
   real(wp) :: c000(1:size(c,4)), c001(1:size(c,4))
   real(wp) :: c010(1:size(c,4)), c011(1:size(c,4))
   real(wp) :: y00, y01, y10, y11
   real(wp) :: y0, y1

   ix1 = int(x(1))
   dx1 = x(1) - real(ix1,wp)

   ix2 = int(x(2))
   dx2 = x(2) - real(ix2,wp)

   ix3 = int(x(3))
   dx3 = x(3) - real(ix3,wp)

   c000 = (1.0_wp-dx1)*c(ix1,ix2,ix3,:) + dx1*c(ix1+1,ix2,ix3,:)
   c001 = (1.0_wp-dx1)*c(ix1,ix2,ix3+1,:) + dx1*c(ix1+1,ix2,ix3+1,:)
   c010 = (1.0_wp-dx1)*c(ix1,ix2+1,ix3,:) + dx1*c(ix1+1,ix2+1,ix3,:)
   c011 = (1.0_wp-dx1)*c(ix1,ix2+1,ix3+1,:) + dx1*c(ix1+1,ix2+1,ix3+1,:)

   y00 = spline1([x(4)], c000, lb, ub)
   y01 = spline1([x(4)], c001, lb, ub)
   y10 = spline1([x(4)], c010, lb, ub)
   y11 = spline1([x(4)], c011, lb, ub)

   y0 = (1.0_wp-dx3)*y00 + dx3*y01
   y1 = (1.0_wp-dx3)*y10 + dx3*y11

   y_interp = (1.0_wp-dx2)*y0 + dx2*y1
end function spline4

function spline4_grid(x, c, lb, ub) result(y_interp)
   real(wp), intent(in) :: x(:)
   real(wp), intent(in) :: c(:,:,:,:)
   real(wp), intent(in) :: lb(:), ub(:)
   real(wp) :: y_interp

   y_interp = spline4(x, c, lb, ub)
end function spline4_grid

!===============================================================================
! 5D Spline Evaluation
!===============================================================================
function spline5(x, c, lb, ub) result(y_interp)
   real(wp), intent(in) :: x(:)
   real(wp), intent(in) :: c(:,:,:,:,:)
   real(wp), intent(in) :: lb(:), ub(:)
   real(wp) :: y_interp

   ! Placeholder — extend in the same pattern as spline4 if needed
   y_interp = 0.0_wp
end function spline5

function spline5_grid(x, c, lb, ub) result(y_interp)
   real(wp), intent(in) :: x(:)
   real(wp), intent(in) :: c(:,:,:,:,:)
   real(wp), intent(in) :: lb(:), ub(:)
   real(wp) :: y_interp

   y_interp = spline5(x, c, lb, ub)
end function spline5_grid

!===============================================================================
! 6D Spline Evaluation (Placeholder)
!===============================================================================
function spline6(x, c, lb, ub) result(y_interp)
   real(wp), intent(in) :: x(:)
   real(wp), intent(in) :: c(:,:,:,:,:,:)
   real(wp), intent(in) :: lb(:), ub(:)
   real(wp) :: y_interp

   ! Placeholder — extend similarly if needed
   y_interp = 0.0_wp
end function spline6

function spline6_grid(x, c, lb, ub) result(y_interp)
   real(wp), intent(in) :: x(:)
   real(wp), intent(in) :: c(:,:,:,:,:,:)
   real(wp), intent(in) :: lb(:), ub(:)
   real(wp) :: y_interp

   y_interp = spline6(x, c, lb, ub)
end function spline6_grid

!===============================================================================
! 7D Spline Evaluation (Placeholder)
!===============================================================================
function spline7(x, c, lb, ub) result(y_interp)
   real(wp), intent(in) :: x(:)
   real(wp), intent(in) :: c(:,:,:,:,:,:,:)
   real(wp), intent(in) :: lb(:), ub(:)
   real(wp) :: y_interp

   ! Placeholder — extend similarly if needed
   y_interp = 0.0_wp
end function spline7

function spline7_grid(x, c, lb, ub) result(y_interp)
   real(wp), intent(in) :: x(:)
   real(wp), intent(in) :: c(:,:,:,:,:,:,:)
   real(wp), intent(in) :: lb(:), ub(:)
   real(wp) :: y_interp

   y_interp = spline7(x, c, lb, ub)
end function spline7_grid

!===============================================================================
! End of Module
!===============================================================================
end module my_spline_mod
