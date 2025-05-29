program test_diagnostics_html
    use info_window
    implicit none
    character(len=100), dimension(15) :: diagnostics
    character(len=256) :: command
    integer :: i, status

    ! Populate diagnostics with test data
    do i = 1, size(diagnostics)
        write(diagnostics(i), '(A, I0)') 'Test diagnostic line ', i
    end do

    ! Call the subroutine to write the HTML file
    call write_diagnostics_to_html("diagnostics.html", diagnostics)

    ! Print a message for testing
    print *, "Diagnostics HTML file created."

    ! Formulate the system command to open the file in Firefox
    command = "firefox diagnostics.html &"

    ! Execute the system call
    call execute_command_line(trim(command), exitstat=status)

    if (status /= 0) then
        print *, "Failed to open the file in Firefox. Please open 'diagnostics.html' manually."
    end if
end program test_diagnostics_html
