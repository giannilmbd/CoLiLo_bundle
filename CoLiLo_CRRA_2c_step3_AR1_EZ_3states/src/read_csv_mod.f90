module read_csv_mod
    use my_kinds_mod, only: wp
    implicit none

contains

    subroutine read_csv(filename, header, data, nlines, rownames)
        use csv_module !JWilliams' code

        implicit none

        character(len=*), intent(in) :: filename

        integer, intent(out), optional :: nlines
        integer :: n_cols, n_rows, n_rh
        real(kind=wp), intent(inout) :: data(:, :)
        character(len=100), intent(out), dimension(size(data, 2)), optional :: header
        character(len=100), intent(out), dimension(size(data, 1)), optional :: rownames
        character(len=100), dimension(:), allocatable :: header_a
        character(len=200) :: fmt
        integer :: i, j, iostat, line, unit
        type(csv_file) :: f
        logical :: status_ok
        ! integer,dimension(:),allocatable :: itypes
        real(kind=wp), dimension(:), allocatable :: tmp_data
        character(len=100), dimension(:), allocatable :: tmp_data_char

        ! character(len=100) ::
        n_rows = size(data, 1)
        n_cols = size(data, 2)
        unit = 10
        n_rh = 0
        if (present(rownames)) then
            n_rh = 1
        end if
        if (present(rownames)) then
            allocate (tmp_data_char(n_rows))
        else
            allocate (tmp_data_char(1))
        end if

        allocate (tmp_data(n_rows))

! read the file
        if (present(header)) then
            call f%read(filename, header_row=1, status_ok=status_ok)
        else
            call f%read(filename, status_ok=status_ok)
        end if

! get the header and type info
        if (present(header)) then
            allocate (header_a(size(data, 2) + n_rh))
            call f%get_header(header_a, status_ok)
            header = header_a(1 + n_rh:)
        else
            allocate (header_a(1))
        end if
! GET ROWNAMES IF EXIST
        if (present(rownames)) then
            call f%get(1, tmp_data_char, status_ok)
            rownames = tmp_data_char
        end if
! call f%variable_types(itypes,status_ok)
        do i = 1 , n_cols
! get some data
            call f%get(i+ n_rh, tmp_data, status_ok)
            if (status_ok .neqv. .true.) then
                print *, '\033[91m PROBLEM IMPORTING DATA '//trim(filename)//'\033[0m'
            end if
! print*,'============READ============'
! print*,tmp_data
! print*,'========================'
            data(:, i) = tmp_data
        end do
! destroy the file
        deallocate (header_a, tmp_data, tmp_data_char)
        call f%destroy()

    end subroutine read_csv

!     ! module save_and_reload_mod
!     !     implicit none

!     !     contains
    subroutine write_csv(filename, header, data, rownames)
        use csv_module !JWilliams' code
        implicit none
        real(kind=wp), intent(in) :: data(:, :)
        character(len=*), intent(in):: filename
        character(len=*), optional, dimension(:) :: header
        character(len=*), optional, dimension(size(data, 1)) :: rownames
        character(len=100), allocatable :: header_ext(:)

        integer :: nlines, n_rows, n_cols, i
        type(csv_file) :: f
        logical :: status_ok

        n_rows = size(data, 1)
        n_cols = size(data, 2)


        if (present(header)) then
            allocate(header_ext(n_cols + 1))
            header_ext(1) = 'rownames'
            header_ext(2:) = header
        end if
        
        ! set optional inputs:
        call f%initialize(verbose=.true.)

        ! open the file
        if (present(rownames)) then
            call f%open(filename, n_cols=n_cols + 1, status_ok=status_ok)
            if (status_ok .neqv. .true.) then
                print *, '\033[91m ERROR IN WRITE_CSV -- FILE '//filename//' NOT FOUND \033[0m'
                stop
            end if
        else
            call f%open(filename, n_cols=n_cols, status_ok=status_ok)
            if (status_ok .neqv. .true.) then
                print *, '\033[91m ERROR IN WRITE_CSV -- FILE '//filename//' NOT FOUND \033[0m'
                stop
            end if
        end if
        ! add header
        if (present(header)) then
        if (present(rownames)) then

            call f%add(header_ext)
        else
            call f%add(header)
        end if
        call f%next_row()
        end if
        deallocate (header_ext)
        ! add some data:
        do i = 1, n_rows
        if (present(rownames)) then
            call f%add(rownames(i), real_fmt='(A)')
        end if
        call f%add(data(i, :), real_fmt='(E20.13)')
        ! print*,'===========WRITE============='
        ! write(*,'(2E20.10)'),data(i,:)
        ! print*,'========================'
        ! call f%add(.true.)
        call f%next_row()
        end do

        ! finished
        call f%close(status_ok)

        !open(unit=10, file=trim(filename), status="replace", action="write")
!     if(present(header))then
!         write(10,'(1x,A20,",")') header
!     endif
! do i = 1, n_rows
!     write(10, '(1x, E20.10, ",")') data(i, 1:n_cols - 1)
!     write(10, '(1x,E20.10)') data(i, n_cols)  ! Last column without a comma
!     ! write(10, *)  ! New line
! end do
! close(10)
    end subroutine write_csv

! ! for debugging (written by Chat GPT)

!     subroutine save_and_reload_array(input_array, discrepancy)
!         real, intent(in) :: input_array(:,:,:,:,:)
!         real, intent(out) :: discrepancy
!         integer, dimension(5) :: array_shape
!         integer :: n_rows, n_cols, i, j
!         real, allocatable :: reshaped_array(:,:), read_array(:,:), temp_array(:,:,:,:,:)

!         ! Get the shape of the input array
!         array_shape = shape(input_array)

!         ! Calculate the size of the reshaped 2D array
!         n_rows = product(array_shape)
!         n_cols = 1

!         ! Reshape the 5D array into a 2D array
!         allocate(reshaped_array(n_rows, n_cols))
!         reshaped_array = reshape(input_array, [n_rows, n_cols])

!         ! Save the 2D array to a CSV file
!         open(unit=10, file="data.csv", status="replace", action="write")
!         do i = 1, n_rows
!             write(10, '(E20.10)') reshaped_array(i, 1)
!         end do
!         close(10)

!         ! Read the 2D array back from the CSV file
!         allocate(read_array(n_rows, n_cols))
!         open(unit=20, file="data.csv", status="old", action="read")
!         do i = 1, n_rows
!             read(20, '(E20.10)') read_array(i, 1)
!         end do
!         close(20)

!         ! Reshape the 2D read array back to the original 5D shape
!         allocate(temp_array(array_shape(1),array_shape(2),array_shape(3),array_shape(4),array_shape(5)))
!         temp_array = reshape(read_array, array_shape)

!         ! Calculate the discrepancy between the original and reloaded arrays
!         discrepancy = maxval(abs(input_array - temp_array))
!     deallocate(temp_array,read_array,reshaped_array)
!     end subroutine save_and_reload_array

!     ! end module save_and_reload_mod

!     ! program main
!     !     use save_and_reload_mod
!     !     implicit none
!     !     real, dimension(2, 2, 2, 2, 2) :: my_array
!     !     real :: discrepancy

!     !     ! Generate some random data for the 5D array
!     !     call random_number(my_array)

!     !     ! Use the subroutine to save and reload the array
!     !     call save_and_reload_array(my_array, discrepancy)

!     !     ! Print the discrepancy
!     !     print *, "Discrepancy:", discrepancy

!     ! end program main

end module read_csv_mod
