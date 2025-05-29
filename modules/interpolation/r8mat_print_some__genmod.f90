        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 18 20:44:15 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE R8MAT_PRINT_SOME__genmod
          INTERFACE 
            RECURSIVE SUBROUTINE R8MAT_PRINT_SOME(M,N,A,ILO,JLO,IHI,JHI,&
     &TITLE)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: M
              REAL(KIND=8) :: A(M,N)
              INTEGER(KIND=4) :: ILO
              INTEGER(KIND=4) :: JLO
              INTEGER(KIND=4) :: IHI
              INTEGER(KIND=4) :: JHI
              CHARACTER(*) :: TITLE
            END SUBROUTINE R8MAT_PRINT_SOME
          END INTERFACE 
        END MODULE R8MAT_PRINT_SOME__genmod
