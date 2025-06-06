# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/blas/crotg.f"
      SUBROUTINE CROTG(CA,CB,C,S)
*     .. Scalar Arguments ..
      COMPLEX CA,CB,S
      REAL C
*     ..
*
*  Purpose
*  =======
*
*  CROTG determines a complex Givens rotation.
*
*  =====================================================================
*
*     .. Local Scalars ..
      COMPLEX ALPHA
      REAL NORM,SCALE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC CABS,CONJG,SQRT
*     ..
      IF (CABS(CA).EQ.0.) THEN
         C = 0.
         S = (1.,0.)
         CA = CB
      ELSE
         SCALE = CABS(CA) + CABS(CB)
         NORM = SCALE*SQRT((CABS(CA/SCALE))**2+ (CABS(CB/SCALE))**2)
         ALPHA = CA/CABS(CA)
         C = CABS(CA)/NORM
         S = ALPHA*CONJG(CB)/NORM
         CA = ALPHA*NORM
      END IF
      RETURN
      END
