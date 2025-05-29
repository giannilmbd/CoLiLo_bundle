        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 18 20:44:15 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE U_POLYNOMIAL_AB__genmod
          INTERFACE 
            RECURSIVE SUBROUTINE U_POLYNOMIAL_AB(A,B,M,N,XAB,V)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: M
              REAL(KIND=8) :: A
              REAL(KIND=8) :: B
              REAL(KIND=8) :: XAB(1:M)
              REAL(KIND=8) :: V(1:M,0:N)
            END SUBROUTINE U_POLYNOMIAL_AB
          END INTERFACE 
        END MODULE U_POLYNOMIAL_AB__genmod
