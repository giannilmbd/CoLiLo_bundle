        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 18 20:44:14 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DAXPY__genmod
          INTERFACE 
            RECURSIVE SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: DA
              REAL(KIND=8) :: DX(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=8) :: DY(*)
              INTEGER(KIND=4) :: INCY
            END SUBROUTINE DAXPY
          END INTERFACE 
        END MODULE DAXPY__genmod
