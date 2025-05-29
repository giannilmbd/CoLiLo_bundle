        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr 18 20:44:15 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SVD_SOLVE__genmod
          INTERFACE 
            RECURSIVE SUBROUTINE SVD_SOLVE(M,N,A,B,X)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: M
              REAL(KIND=8) :: A(M,N)
              REAL(KIND=8) :: B(M)
              REAL(KIND=8) :: X(N)
            END SUBROUTINE SVD_SOLVE
          END INTERFACE 
        END MODULE SVD_SOLVE__genmod
