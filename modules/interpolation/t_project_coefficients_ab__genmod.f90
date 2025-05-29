        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 18 20:44:15 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE T_PROJECT_COEFFICIENTS_AB__genmod
          INTERFACE 
            RECURSIVE SUBROUTINE T_PROJECT_COEFFICIENTS_AB(N,F,A,B,C)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: F
              EXTERNAL F
              REAL(KIND=8) :: A
              REAL(KIND=8) :: B
              REAL(KIND=8) :: C(0:N)
            END SUBROUTINE T_PROJECT_COEFFICIENTS_AB
          END INTERFACE 
        END MODULE T_PROJECT_COEFFICIENTS_AB__genmod
