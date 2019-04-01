!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Constants, Data, and functions for OSB package in sequential environment   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  module OSBW
      use OSB
      implicit none
      integer, parameter :: MAXPOLYNUM=20
      double precision :: polyCoeff(MAXPOLYNUM+1)
      integer :: polyNum

      contains   ! Function defined for OSB/OSBW

        include 'osbw.qmr.h'        ! QMR method
        include 'osbw.qmrdx.h'      !
        include 'osbw.qmrcx.h'      !   

!        include 'osbw.hij.h'        

        include 'osbw.lanczos.h'    ! Lanczos algorithm, only for real/Hermitian
        include 'osbw.polylan.h'    ! Polynomial Lanczos, only for real/Hermitian
        include 'osbw.pist.h'       ! PIST algorithm, real/complex
        include 'osbw.pistconv.h'
 
             ! CRP calculation
        include 'osbw.crp.h'        ! CRP calculation,

             ! testing subroutines
        include 'osbw.progQMR.h'
        include 'osbw.progEig.h' 
!        include 'osbw.progLan.h' 

   end module OSBW
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
