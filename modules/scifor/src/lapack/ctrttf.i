# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/lapack/ctrttf.f"
      SUBROUTINE CTRTTF( TRANSR, UPLO, N, A, LDA, ARF, INFO )
*
*  -- LAPACK routine (version 3.3.1)                                    --
*
*  -- Contributed by Fred Gustavson of the IBM Watson Research Center --
*  -- April 2011                                                      --
*
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          TRANSR, UPLO
      INTEGER            INFO, N, LDA
*     ..
*     .. Array Arguments ..
      COMPLEX            A( 0: LDA-1, 0: * ), ARF( 0: * )
*     ..
*
*  Purpose
*  =======
*
*  CTRTTF copies a triangular matrix A from standard full format (TR)
*  to rectangular full packed format (TF) .
*
*  Arguments
*  =========
*
*  TRANSR   (input) CHARACTER*1
*          = 'N':  ARF in Normal mode is wanted;
*          = 'C':  ARF in Conjugate Transpose mode is wanted;
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  A is upper triangular;
*          = 'L':  A is lower triangular.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input) COMPLEX array, dimension ( LDA, N )
*          On entry, the triangular matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of the array A contains
*          the upper triangular matrix, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of the array A contains
*          the lower triangular matrix, and the strictly upper
*          triangular part of A is not referenced.
*
*  LDA     (input) INTEGER
*          The leading dimension of the matrix A.  LDA >= max(1,N).
*
*  ARF     (output) COMPLEX*16 array, dimension ( N*(N+1)/2 ),
*          On exit, the upper or lower triangular matrix A stored in
*          RFP format. For a further discussion see Notes below.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  We first consider Standard Packed Format when N is even.
*  We give an example where N = 6.
*
*      AP is Upper             AP is Lower
*
*   00 01 02 03 04 05       00
*      11 12 13 14 15       10 11
*         22 23 24 25       20 21 22
*            33 34 35       30 31 32 33
*               44 45       40 41 42 43 44
*                  55       50 51 52 53 54 55
*
*
*  Let TRANSR = 'N'. RFP holds AP as follows:
*  For UPLO = 'U' the upper trapezoid A(0:5,0:2) consists of the last
*  three columns of AP upper. The lower triangle A(4:6,0:2) consists of
*  conjugate-transpose of the first three columns of AP upper.
*  For UPLO = 'L' the lower trapezoid A(1:6,0:2) consists of the first
*  three columns of AP lower. The upper triangle A(0:2,0:2) consists of
*  conjugate-transpose of the last three columns of AP lower.
*  To denote conjugate we place -- above the element. This covers the
*  case N even and TRANSR = 'N'.
*
*         RFP A                   RFP A
*
*                                -- -- --
*        03 04 05                33 43 53
*                                   -- --
*        13 14 15                00 44 54
*                                      --
*        23 24 25                10 11 55
*
*        33 34 35                20 21 22
*        --
*        00 44 45                30 31 32
*        -- --
*        01 11 55                40 41 42
*        -- -- --
*        02 12 22                50 51 52
*
*  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate-
*  transpose of RFP A above. One therefore gets:
*
*
*           RFP A                   RFP A
*
*     -- -- -- --                -- -- -- -- -- --
*     03 13 23 33 00 01 02    33 00 10 20 30 40 50
*     -- -- -- -- --                -- -- -- -- --
*     04 14 24 34 44 11 12    43 44 11 21 31 41 51
*     -- -- -- -- -- --                -- -- -- --
*     05 15 25 35 45 55 22    53 54 55 22 32 42 52
*
*
*  We next  consider Standard Packed Format when N is odd.
*  We give an example where N = 5.
*
*     AP is Upper                 AP is Lower
*
*   00 01 02 03 04              00
*      11 12 13 14              10 11
*         22 23 24              20 21 22
*            33 34              30 31 32 33
*               44              40 41 42 43 44
*
*
*  Let TRANSR = 'N'. RFP holds AP as follows:
*  For UPLO = 'U' the upper trapezoid A(0:4,0:2) consists of the last
*  three columns of AP upper. The lower triangle A(3:4,0:1) consists of
*  conjugate-transpose of the first two   columns of AP upper.
*  For UPLO = 'L' the lower trapezoid A(0:4,0:2) consists of the first
*  three columns of AP lower. The upper triangle A(0:1,1:2) consists of
*  conjugate-transpose of the last two   columns of AP lower.
*  To denote conjugate we place -- above the element. This covers the
*  case N odd  and TRANSR = 'N'.
*
*         RFP A                   RFP A
*
*                                   -- --
*        02 03 04                00 33 43
*                                      --
*        12 13 14                10 11 44
*
*        22 23 24                20 21 22
*        --
*        00 33 34                30 31 32
*        -- --
*        01 11 44                40 41 42
*
*  Now let TRANSR = 'C'. RFP A in both UPLO cases is just the conjugate-
*  transpose of RFP A above. One therefore gets:
*
*
*           RFP A                   RFP A
*
*     -- -- --                   -- -- -- -- -- --
*     02 12 22 00 01             00 10 20 30 40 50
*     -- -- -- --                   -- -- -- -- --
*     03 13 23 33 11             33 11 21 31 41 51
*     -- -- -- -- --                   -- -- -- --
*     04 14 24 34 44             43 44 22 32 42 52
*
*  =====================================================================
*
*     .. Parameters ..
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, NISODD, NORMALTRANSR
      INTEGER            I, IJ, J, K, L, N1, N2, NT, NX2, NP1X2
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX, MOD
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NORMALTRANSR = LSAME( TRANSR, 'N' )
      LOWER = LSAME( UPLO, 'L' )
      IF( .NOT.NORMALTRANSR .AND. .NOT.LSAME( TRANSR, 'C' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTRTTF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.1 ) THEN
         IF( N.EQ.1 ) THEN
            IF( NORMALTRANSR ) THEN
               ARF( 0 ) = A( 0, 0 )
            ELSE
               ARF( 0 ) = CONJG( A( 0, 0 ) )
            END IF
         END IF
         RETURN
      END IF
*
*     Size of array ARF(1:2,0:nt-1)
*
      NT = N*( N+1 ) / 2
*
*     set N1 and N2 depending on LOWER: for N even N1=N2=K
*
      IF( LOWER ) THEN
         N2 = N / 2
         N1 = N - N2
      ELSE
         N1 = N / 2
         N2 = N - N1
      END IF
*
*     If N is odd, set NISODD = .TRUE., LDA=N+1 and A is (N+1)--by--K2.
*     If N is even, set K = N/2 and NISODD = .FALSE., LDA=N and A is
*     N--by--(N+1)/2.
*
      IF( MOD( N, 2 ).EQ.0 ) THEN
         K = N / 2
         NISODD = .FALSE.
         IF( .NOT.LOWER )
     $      NP1X2 = N + N + 2
      ELSE
         NISODD = .TRUE.
         IF( .NOT.LOWER )
     $      NX2 = N + N
      END IF
*
      IF( NISODD ) THEN
*
*        N is odd
*
         IF( NORMALTRANSR ) THEN
*
*           N is odd and TRANSR = 'N'
*
            IF( LOWER ) THEN
*
*             SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) )
*             T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0)
*             T1 -> a(0), T2 -> a(n), S -> a(n1); lda=n
*
               IJ = 0
               DO J = 0, N2
                  DO I = N1, N2 + J
                     ARF( IJ ) = CONJG( A( N2+J, I ) )
                     IJ = IJ + 1
                  END DO
                  DO I = J, N - 1
                     ARF( IJ ) = A( I, J )
                     IJ = IJ + 1
                  END DO
               END DO
*
            ELSE
*
*             SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
*             T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
*             T1 -> a(n2), T2 -> a(n1), S -> a(0); lda=n
*
               IJ = NT - N
               DO J = N - 1, N1, -1
                  DO I = 0, J
                     ARF( IJ ) = A( I, J )
                     IJ = IJ + 1
                  END DO
                  DO L = J - N1, N1 - 1
                     ARF( IJ ) = CONJG( A( J-N1, L ) )
                     IJ = IJ + 1
                  END DO
                  IJ = IJ - NX2
               END DO
*
            END IF
*
         ELSE
*
*           N is odd and TRANSR = 'C'
*
            IF( LOWER ) THEN
*
*              SRPA for LOWER, TRANSPOSE and N is odd
*              T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1)
*              T1 -> A(0+0) , T2 -> A(1+0) , S -> A(0+n1*n1); lda=n1
*
               IJ = 0
               DO J = 0, N2 - 1
                  DO I = 0, J
                     ARF( IJ ) = CONJG( A( J, I ) )
                     IJ = IJ + 1
                  END DO
                  DO I = N1 + J, N - 1
                     ARF( IJ ) = A( I, N1+J )
                     IJ = IJ + 1
                  END DO
               END DO
               DO J = N2, N - 1
                  DO I = 0, N1 - 1
                     ARF( IJ ) = CONJG( A( J, I ) )
                     IJ = IJ + 1
                  END DO
               END DO
*
            ELSE
*
*              SRPA for UPPER, TRANSPOSE and N is odd
*              T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0)
*              T1 -> A(n2*n2), T2 -> A(n1*n2), S -> A(0); lda=n2
*
               IJ = 0
               DO J = 0, N1
                  DO I = N1, N - 1
                     ARF( IJ ) = CONJG( A( J, I ) )
                     IJ = IJ + 1
                  END DO
               END DO
               DO J = 0, N1 - 1
                  DO I = 0, J
                     ARF( IJ ) = A( I, J )
                     IJ = IJ + 1
                  END DO
                  DO L = N2 + J, N - 1
                     ARF( IJ ) = CONJG( A( N2+J, L ) )
                     IJ = IJ + 1
                  END DO
               END DO
*
            END IF
*
         END IF
*
      ELSE
*
*        N is even
*
         IF( NORMALTRANSR ) THEN
*
*           N is even and TRANSR = 'N'
*
            IF( LOWER ) THEN
*
*              SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) )
*              T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0)
*              T1 -> a(1), T2 -> a(0), S -> a(k+1); lda=n+1
*
               IJ = 0
               DO J = 0, K - 1
                  DO I = K, K + J
                     ARF( IJ ) = CONJG( A( K+J, I ) )
                     IJ = IJ + 1
                  END DO
                  DO I = J, N - 1
                     ARF( IJ ) = A( I, J )
                     IJ = IJ + 1
                  END DO
               END DO
*
            ELSE
*
*              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
*              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
*              T1 -> a(k+1), T2 -> a(k), S -> a(0); lda=n+1
*
               IJ = NT - N - 1
               DO J = N - 1, K, -1
                  DO I = 0, J
                     ARF( IJ ) = A( I, J )
                     IJ = IJ + 1
                  END DO
                  DO L = J - K, K - 1
                     ARF( IJ ) = CONJG( A( J-K, L ) )
                     IJ = IJ + 1
                  END DO
                  IJ = IJ - NP1X2
               END DO
*
            END IF
*
         ELSE
*
*           N is even and TRANSR = 'C'
*
            IF( LOWER ) THEN
*
*              SRPA for LOWER, TRANSPOSE and N is even (see paper, A=B)
*              T1 -> A(0,1) , T2 -> A(0,0) , S -> A(0,k+1) :
*              T1 -> A(0+k) , T2 -> A(0+0) , S -> A(0+k*(k+1)); lda=k
*
               IJ = 0
               J = K
               DO I = K, N - 1
                  ARF( IJ ) = A( I, J )
                  IJ = IJ + 1
               END DO
               DO J = 0, K - 2
                  DO I = 0, J
                     ARF( IJ ) = CONJG( A( J, I ) )
                     IJ = IJ + 1
                  END DO
                  DO I = K + 1 + J, N - 1
                     ARF( IJ ) = A( I, K+1+J )
                     IJ = IJ + 1
                  END DO
               END DO
               DO J = K - 1, N - 1
                  DO I = 0, K - 1
                     ARF( IJ ) = CONJG( A( J, I ) )
                     IJ = IJ + 1
                  END DO
               END DO
*
            ELSE
*
*              SRPA for UPPER, TRANSPOSE and N is even (see paper, A=B)
*              T1 -> A(0,k+1) , T2 -> A(0,k) , S -> A(0,0)
*              T1 -> A(0+k*(k+1)) , T2 -> A(0+k*k) , S -> A(0+0)); lda=k
*
               IJ = 0
               DO J = 0, K
                  DO I = K, N - 1
                     ARF( IJ ) = CONJG( A( J, I ) )
                     IJ = IJ + 1
                  END DO
               END DO
               DO J = 0, K - 2
                  DO I = 0, J
                     ARF( IJ ) = A( I, J )
                     IJ = IJ + 1
                  END DO
                  DO L = K + 1 + J, N - 1
                     ARF( IJ ) = CONJG( A( K+1+J, L ) )
                     IJ = IJ + 1
                  END DO
               END DO
*
*              Note that here J = K-1
*
               DO I = 0, J
                  ARF( IJ ) = A( I, J )
                  IJ = IJ + 1
               END DO
*
            END IF
*
         END IF
*
      END IF
*
      RETURN
*
*     End of CTRTTF
*
      END
