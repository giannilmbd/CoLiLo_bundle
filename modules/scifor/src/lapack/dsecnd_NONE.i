# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/dsecnd_NONE.f"
      DOUBLE PRECISION FUNCTION DSECND( )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     February 2007
*
*  Purpose
*  =======
*
*  DSECND returns nothing instead of returning the user time for a process in seconds.
*  If you are using that routine, it means that neither EXTERNAL ETIME,
*  EXTERNAL ETIME_, INTERNAL ETIME, INTERNAL CPU_TIME is available  on
*  your machine.
*
* =====================================================================
*
      DSECND = 0.0D+0
      RETURN
*
*     End of DSECND
*
      END
