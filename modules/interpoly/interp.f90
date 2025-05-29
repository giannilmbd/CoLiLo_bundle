program interp
    use fgsl
    implicit none
    integer, parameter :: order = 4
    integer :: i
    integer, parameter :: n = 5
    real(fgsl_double) :: xi(n), yi(n), x(100), y(100)
    type(fgsl_interp_type_ptr) :: interp_type_ptr
    type(fgsl_interp) :: interp
    type(fgsl_interp_accel) :: accel

    ! set up x and y arrays for interpolation
    xi = [0.0_fgsl_double, 0.2_fgsl_double, 0.4_fgsl_double, 0.6_fgsl_double, 0.8_fgsl_double]
    yi = [0.0_fgsl_double, 0.79465447229_fgsl_double, 1.30540728979_fgsl_double, &
          1.57079632679_fgsl_double, 1.32460908925_fgsl_double]

    ! set up x array for evaluation
    x = [(0.01_fgsl_double * real(i, kind=fgsl_double)), i = 0, 99]

    ! create interpolation type object
    interp_type_ptr => fgsl_interp_polynomial

    ! create interpolation object and accelerator
    interp = fgsl_interp_alloc(interp_type_ptr, n)
    accel = fgsl_interp_accel_alloc()

    ! initialize interpolation object
    call fgsl_interp_init(interp, xi, yi, n)

    ! interpolate y values for x values
    do i = 1, 100
        y(i) = fgsl_interp_eval(interp, x(i), accel)
    end do

    ! print out interpolated values
    do i = 1, 100
        write(*, '(F6.2, 2X, F10.6)') x(i), y(i)
    end do

    ! free allocated memory
    call fgsl_interp_free(interp)
    call fgsl_interp_accel_free(accel)
end program interp

  
  !ifort interp.f90 -I/usr/local/include/fgsl -L/usr/local/lib/ -lfgsl   -o interp.out