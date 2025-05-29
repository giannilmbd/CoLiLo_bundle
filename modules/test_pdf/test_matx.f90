program test_matx
    use simulated_annealing_module, only: print_matrix
    implicit none
   
real*8 ,dimension(5,5) ::matx
integer :: cnt,sz=5
character*20 :: txt    

matx(:,:)=1.01
call print_matrix(sz,sz,matx,'text')

end program test_matx