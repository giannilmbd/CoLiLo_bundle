# 1 "/home/gianni/Dropbox/projects/ANNAGIANCARLOGIANNI/CoLiLoGlobal/modules/mylib/rantest.f"
      program rantest
C
C  illustrates use of function RAN3
C  computes the average
C
C  CWJ SDSU 2/20/2005
C
C  functions called: ran3 (found in ranlib.f)
C
      implicit none

      real ran3
      integer iseed
      real avg,avg2
      integer nsample
      real x
      integer i

      print*,' please enter an integer '
      read*,iseed
      if(iseed.gt.0)iseed=-iseed
      print*,' How many samples ?'
      read*,nsample

      avg = 0.0

      do i = 1,nsample
        x=ran3(iseed)
        print*,x
        avg=avg+x
        avg2 = avg2+x*x
      enddo
 
      avg = avg/nsample

      avg2 = avg2/nsample-(avg)**2

      write(6,101)avg,sqrt(avg2/nsample)
  101 format(' avg = ',f5.3,' +/- ',f5.3)


      end

