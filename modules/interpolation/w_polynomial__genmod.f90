        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 18 20:44:15 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE W_POLYNOMIAL__genmod
          INTERFACE 
            RECURSIVE SUBROUTINE W_POLYNOMIAL(M,N,X,V)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: M
              REAL(KIND=8) :: X(M)
              REAL(KIND=8) :: V(M,0:N)
            END SUBROUTINE W_POLYNOMIAL
          END INTERFACE 
        END MODULE W_POLYNOMIAL__genmod
