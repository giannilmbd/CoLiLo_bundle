        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 18 20:44:15 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DSVDC__genmod
          INTERFACE 
            RECURSIVE SUBROUTINE DSVDC(A,LDA,M,N,S,E,U,LDU,V,LDV,WORK,  &
     &JOB,INFO)
              INTEGER(KIND=4) :: LDV
              INTEGER(KIND=4) :: LDU
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: LDA
              REAL(KIND=8) :: A(LDA,N)
              REAL(KIND=8) :: S(*)
              REAL(KIND=8) :: E(*)
              REAL(KIND=8) :: U(LDU,M)
              REAL(KIND=8) :: V(LDV,N)
              REAL(KIND=8) :: WORK(M)
              INTEGER(KIND=4) :: JOB
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DSVDC
          END INTERFACE 
        END MODULE DSVDC__genmod
