        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 18 20:44:15 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE R8MAT_PRINT__genmod
          INTERFACE 
            RECURSIVE SUBROUTINE R8MAT_PRINT(M,N,A,TITLE)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: M
              REAL(KIND=8) :: A(M,N)
              CHARACTER(*) :: TITLE
            END SUBROUTINE R8MAT_PRINT
          END INTERFACE 
        END MODULE R8MAT_PRINT__genmod
