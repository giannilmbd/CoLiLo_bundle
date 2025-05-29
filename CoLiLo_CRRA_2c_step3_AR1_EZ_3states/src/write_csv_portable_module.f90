module write_csv_portable_module
   
   implicit none

contains
   ! This subroutine writes a 2D array of real numbers to a CSV file.
   ! It uses the Fortran intrinsic I/O capabilities to ensure compatibility
   ! across different platforms and compilers.
   ! The data is written in scientific notation with 13 decimal places.

   ! Note: The file is opened with 'replace' status, which means it will
   ! overwrite any existing file with the same name.
subroutine write_csv_portable(filename, data)
   use my_kinds_mod
   implicit none
   character(len=*), intent(in) :: filename
   real(kind=wp), intent(in)    :: data(:, :)
   integer :: i, j, n_rows, n_cols
   integer :: unit

   n_rows = size(data, 1)
   n_cols = size(data, 2)

   ! Try to open the file
   open(newunit=unit, file=filename, status='replace', action='write', iostat=i)
   if (i /= 0) then
      print *, 'ERROR: Unable to open file for writing:', trim(filename)
      stop 1
   end if

   ! Write data row by row
   do i = 1, n_rows
      do j = 1, n_cols
         write(unit, '(E20.13)', advance='no') data(i, j)
         if (j < n_cols) write(unit, '(A)', advance='no') ','
      end do
      write(unit, *)
   end do

   close(unit)
end subroutine write_csv_portable
end module write_csv_portable_module