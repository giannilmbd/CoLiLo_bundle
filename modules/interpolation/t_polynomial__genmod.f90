        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 18 20:44:15 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE T_POLYNOMIAL__genmod
          INTERFACE 
            RECURSIVE SUBROUTINE T_POLYNOMIAL(M,N,X,V)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: M
              REAL(KIND=8) :: X(1:M)
              REAL(KIND=8) :: V(1:M,0:N)
            END SUBROUTINE T_POLYNOMIAL
          END INTERFACE 
        END MODULE T_POLYNOMIAL__genmod
