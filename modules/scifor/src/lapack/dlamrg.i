# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/dlamrg.f"
      SUBROUTINE DLAMRG( N1, N2, A, DTRD1, DTRD2, INDEX )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            DTRD1, DTRD2, N1, N2
*     ..
*     .. Array Arguments ..
      INTEGER            INDEX( * )
      DOUBLE PRECISION   A( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAMRG will create a permutation list which will merge the elements
*  of A (which is composed of two independently sorted sets) into a
*  single set which is sorted in ascending order.
*
*  Arguments
*  =========
*
*  N1     (input) INTEGER
*  N2     (input) INTEGER
*         These arguements contain the respective lengths of the two
*         sorted lists to be merged.
*
*  A      (input) DOUBLE PRECISION array, dimension (N1+N2)
*         The first N1 elements of A contain a list of numbers which
*         are sorted in either ascending or descending order.  Likewise
*         for the final N2 elements.
*
*  DTRD1  (input) INTEGER
*  DTRD2  (input) INTEGER
*         These are the strides to be taken through the array A.
*         Allowable strides are 1 and -1.  They indicate whether a
*         subset of A is sorted in ascending (DTRDx = 1) or descending
*         (DTRDx = -1) order.
*
*  INDEX  (output) INTEGER array, dimension (N1+N2)
*         On exit this array will contain a permutation such that
*         if B( I ) = A( INDEX( I ) ) for I=1,N1+N2, then B will be
*         sorted in ascending order.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IND1, IND2, N1SV, N2SV
*     ..
*     .. Executable Statements ..
*
      N1SV = N1
      N2SV = N2
      IF( DTRD1.GT.0 ) THEN
         IND1 = 1
      ELSE
         IND1 = N1
      END IF
      IF( DTRD2.GT.0 ) THEN
         IND2 = 1 + N1
      ELSE
         IND2 = N1 + N2
      END IF
      I = 1
*     while ( (N1SV > 0) & (N2SV > 0) )
   10 CONTINUE
      IF( N1SV.GT.0 .AND. N2SV.GT.0 ) THEN
         IF( A( IND1 ).LE.A( IND2 ) ) THEN
            INDEX( I ) = IND1
            I = I + 1
            IND1 = IND1 + DTRD1
            N1SV = N1SV - 1
         ELSE
            INDEX( I ) = IND2
            I = I + 1
            IND2 = IND2 + DTRD2
            N2SV = N2SV - 1
         END IF
         GO TO 10
      END IF
*     end while
      IF( N1SV.EQ.0 ) THEN
         DO 20 N1SV = 1, N2SV
            INDEX( I ) = IND2
            I = I + 1
            IND2 = IND2 + DTRD2
   20    CONTINUE
      ELSE
*     N2SV .EQ. 0
         DO 30 N2SV = 1, N1SV
            INDEX( I ) = IND1
            I = I + 1
            IND1 = IND1 + DTRD1
   30    CONTINUE
      END IF
*
      RETURN
*
*     End of DLAMRG
*
      END
