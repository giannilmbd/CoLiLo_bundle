# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/arpack/util/second.f"
      SUBROUTINE ARSCND( T )
*
      REAL       T
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     July 26, 1991
*
*  Purpose
*  =======
*
*  SECOND returns the user time for a process in seconds.
*  This version gets the time from the system function ETIME.
*
*     .. Local Scalars ..
      REAL               T1
*     ..
*     .. Local Arrays ..
      REAL               TARRAY( 2 )
*     ..
*     .. External Functions ..
      REAL               ETIME
*      INTRINSIC          ETIME
*     ..
*     .. Executable Statements ..
*

      T1 = ETIME( TARRAY )
      T  = TARRAY( 1 )

      RETURN
*
*     End of ARSCND
*
      END
