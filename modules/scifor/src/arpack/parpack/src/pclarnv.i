# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/scifor/src/arpack/parpack/src/pclarnv.f"
c\BeginDoc
c
c\Name: pclarnv
c
c Message Passing Layer: MPI
c
c\Description:
c
c  Parallel Version of ARPACK utility routine clarnv
c
c  PCLARNV returns a vector of n (nloc) random Complex numbers from a uniform or
c  normal distribution. It is assumed that X is distributed across a 1-D array
c  of processors ( nprocs < 1000 )
c
c\Arguments
c  COMM    MPI Communicator for the processor grid
c
c  IDIST   (input) INTEGER
c          Specifies the distribution of the random numbers:
c          = 1:  uniform (0,1)
c          = 2:  uniform (-1,1)
c          = 3:  normal (0,1)
c
c  ISEED   (input/output) INTEGER array, dimension (4)
c          On entry, the seed of the random number generator; the array
c          elements must be between 0 and 4095, and ISEED(4) must be
c          odd.
c          On exit, the seed is updated.
c
c  N       (input) INTEGER
c          The number of random numbers to be generated.
c
c  X       (output) Complex array, dimension (N)
c          The generated random numbers.
c
c\Author: Kristi Maschhoff
c
c\Details
c
c  Simple parallel version of LAPACK auxiliary routine clarnv
c  for X distributed across a 1-D array of processors.
c  This routine calls the auxiliary routine CLARNV to generate random
c  Complex numbers from a uniform or normal distribution. Output is consistent
c  with serial version.
c
c\SCCS Information:
c FILE: larnv.F   SID: 1.3   DATE OF SID: 04/17/99
c
c-----------------------------------------------------------------------
c
      subroutine pclarnv( comm, idist, iseed, n, x )
c
      integer   comm
c     ..
c     .. Scalar Arguments ..
      integer			idist, n
c     ..
c     .. Array Arguments ..
      integer			iseed( 4 )
      Complex
     &                          x( * )
c     ..
c     .. External Subroutines ..
      external			clarnv
c     ..
c     .. Executable Statements ..
c
      call clarnv ( idist, iseed, n, x )
c
      return
      end
