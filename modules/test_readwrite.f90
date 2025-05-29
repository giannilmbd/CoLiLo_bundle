program csv_data_test
    implicit none
    integer, parameter :: n1 = 2, n2 = 2, n3 = 2, n4 = 2, n5 = 2
    real, dimension(n1, n2, n3, n4, n5) :: original_data, final_data
    real, dimension(n1 * n2 * n3 * n4 * n5, 5) :: reshaped_data, read_data
    integer :: i, j, k, l, m, counter
    real :: discrepancy

    ! Generate some data for the 5D array
    do i = 1, n1
        do j = 1, n2
            do k = 1, n3
                do l = 1, n4
                    do m = 1, n5
                        original_data(i, j, k, l, m) = real(i * j * k * l * m)
                    end do
                end do
            end do
        end do
    end do

    ! Reshape the 5D array to a 2D array
    counter = 1
    do i = 1, n1
        do j = 1, n2
            do k = 1, n3
                do l = 1, n4
                    do m = 1, n5
                        reshaped_data(counter, :) = [real(i), real(j), real(k), real(l), real(m), original_data(i, j, k, l, m)]
                        counter = counter + 1
                    end do
                end do
            end do
        end do
    end do

    ! Write the 2D array to a CSV file
    open(10, file='data.csv', status='replace')
    do i = 1, size(reshaped_data, 1)
        write(10, '(5(F0,",",F0),F0)') (reshaped_data(i, j), j = 1, size(reshaped_data, 2))
    end do
    close(10)

    ! Read the 2D array from the CSV file
    open(20, file='data.csv', status='old')
    do i = 1, size(read_data, 1)
        read(20, *, iostat=read_status) (read_data(i, j), j = 1, size(read_data, 2))
        if (read_status /= 0) exit
    end do
    close(20)

    ! Reshape the read 2D array back to a 5D array
    do i = 1, size(read_data, 1)
        final_data(int(read_data(i, 1)), int(read_data(i, 2)), int(read_data(i, 3)), int(read_data(i, 4)), int(read_data(i, 5))) = read_data(i, 6)
    end do

    ! Compute the discrepancy between the original and final 5D arrays
    discrepancy = maxval(abs(final_data - original_data))

    ! Print the discrepancy
    print *, "Discrepancy =", discrepancy

end program csv_data_test
