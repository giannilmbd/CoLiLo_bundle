!     Last change:  GL   16 Feb 2006   10:15 am
! #include "fintrf.h" "matrix.h" "mex.h" "tmwtypes.h"
SUBROUTINE MEXFUNCTION(nlhs,plhs,nrhs,prhs)

  use mexf90                    ! API function definitions

  implicit none
  integer, intent(in) :: nlhs, nrhs
  integer, intent(in), dimension(*) :: prhs
  integer, intent(out), dimension(*) :: plhs
EXTERNAL kron_d  ! the code that want to translate

! LOCAL VARIABLES
 !! used within mexfunction only
! INTEGER, DIMENSION(nrhs) :: rp  ! will be pointers to RHS arguments / there are nrhs RHS arguments so need a nrhs vector
INTEGER :: indxi, nolhs
!double precision, ALLOCATABLE, DIMENSION(:,:) :: a,b,q,z
INTEGER, DIMENSION(nrhs) :: m,n ! rows and cols
! INTEGER, ALLOCATABLE, DIMENSION(:) :: outm, outn, lp !  dimensions and pointers to LHS arguments
integer,pointer :: A,B,AB,periods



nullify(A,B,AB)

! INTEGER :: outp, status!------------------------------------------------
! allocate pointers and size of output arguments
nolhs = max(1,nlhs) ! guarantee that the LHS is at least one object
! ALLOCATE(lp(nolhs),STAT=status) ! left pointer is what just given line above
! !ALLOCATE(lpi(nolhs),STAT=status)
! ALLOCATE(outm(nolhs),STAT=status)
! ALLOCATE(outn(nolhs),STAT=status)
!------------------------------------------------
if (nrhs<4) then  !
  !     call MEXERRMSGTXT('Wrong number of inputs: need 2 matrices')
end if

! CHECK DIMENSION OF Right ARGUMENTS
DO indxi=1,nrhs ! loop through the RHS arguments
      m(indxi) = MXGETM(prhs(indxi)) ! get rows of each of them
      n(indxi) = MXGETN(prhs(indxi)) ! get cols of each of them
END DO

! Read input arguments
! ASSIGN POINTERS TO THE VARIOUS PARAMETERS
! DO indxi=1,nrhs
! 	rp(indxi) = MXGETPR(prhs(indxi)) ! get the memory position of each RHS argument
! END DO



!! Create matrices for return arguments


!ALLOCATE(a(m(1),m(1)))
plhs(1) = MXCREATEDOUBLEMATRIX(m(3),1,0)! size of the output matrix (first output)
!!second output


!





! define pointers to output arguments
! DO indxi =1,nolhs
! 	lp(indxi) = MXGETPR(plhs(indxi))
! !        lpi(indxi) = MXGETPi(plhs(indxi))
! END DO

A => mxGetPr(prhs(1))
B => mxGetPr(prhs(2))
periods=>mxGetPr(prhs(3))

!   lp(1) = mxCreateDoubleMatrix(m,n,0)
  AB=> mxGetPr(plhs(1))

! External business

! CALL kron_d(%VAL(rp(1)),%VAL(rp(2)),m(1),n(1),m(2),n(2),%VAL(lp(1)))
CALL mat_gennorm(A,B,mxGetPr,m(1),n(1),m(3),AB)
!call MXCOPYCOMPLEX16TOPTR(a,lp(1),lpi(1),outn(1)*outm(1))
!call MXCOPYCOMPLEX16TOPTR(b,lp(2),lpi(2),outn(2)*outm(2))
!call MXCOPYCOMPLEX16TOPTR(q,lp(3),lpi(3),outn(3)*outm(3))
!call MXCOPYCOMPLEX16TOPTR(z,lp(4),lpi(4),outn(4)*outm(4))

!DEALLOCATE(a,b,q,z)
! DEALLOCATE(lp,outm,outn);
nullify(A,B,AB)
END SUBROUTINE