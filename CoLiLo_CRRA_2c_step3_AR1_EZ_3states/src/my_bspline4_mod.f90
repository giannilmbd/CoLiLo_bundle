module my_bspline4_mod
  use bspline_oo_module      ! from bspline-fortran
  use bspline_kinds_module, only: wp, ip
  implicit none

  private
  public :: Bspline4, bspline4_build, bspline4_eval

  type :: Bspline4
    private
    type(bspline_4d) :: s4      ! the “native” spline object
  end type Bspline4

contains

  !------------------------------------------------------------
  !> Build a 4-D B-spline from a regular grid + 4-D array of
  !!  function values.  k is the common order in each dimension.
  !!  extrap = .true. allows linear extrapolation at the edges.
  subroutine bspline4_build(grid1,grid2,grid3,grid4,values,k,extrap, spline, iflag)
    real(wp), intent(in)   :: grid1(:), grid2(:), grid3(:), grid4(:)
    real(wp), intent(in)   :: values(size(grid1),size(grid2), &
                                        size(grid3),size(grid4))
    integer(ip), intent(in):: k
    logical,   intent(in), optional :: extrap
    type(Bspline4), intent(out) :: spline
    integer,   intent(out)      :: iflag

    ! call the order-k constructor in each dim
    call spline%s4%initialize( &
         grid1, grid2, grid3, grid4, &
         values, &
         k, k, k, k, &
         iflag, &
         merge(extrap, .false., present(extrap)) )
    if (iflag /= 0) then
      write(*,*) "bspline4_build failed, flag=", iflag
    end if
  end subroutine

  !------------------------------------------------------------
  !> Evaluate the 4-D spline at the point x = x_states(:).  
  !!  Returns the interpolated (or extrapolated) value.
  function bspline4_eval(spline, x_states, iflag) result(fx)
    type(Bspline4), intent(inout)    :: spline
    real(wp),      intent(in)     :: x_states(4)
    integer,       intent(out)    :: iflag
    real(wp)                       :: fx
    integer(ip) :: idx1, idx2, idx3, idx4

    call spline%s4%evaluate( &
      x_states(1), x_states(2), x_states(3), x_states(4), &
      idx1, idx2, idx3, idx4, &
      fx, iflag )

    if (iflag /= 0) then
      write(*,*) "bspline4_eval: evaluation failed, flag=", iflag
    end if
  end function

end module
