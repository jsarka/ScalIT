!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Constants, Data, and functions for OSB package in sequential environment   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  module POSBP
      use POSB
      implicit none
      integer, parameter :: MAXPOLYNUM=20

      double precision :: polyCoeff(MAXPOLYNUM+1)

      integer :: polyNum

      integer :: OSBType

      double precision, allocatable :: SQ_APR(:)

      contains   ! Function defined for OSB/OSBW

        include 'mosbp.qmr.h'        ! QMR method
        include 'mosbp.qmrdx.h'      !
        include 'mosbp.qmrcx.h'      !

        include 'mosbp.lan.h'    ! Lanczos algorithm, only for real/Hermitian
        include 'mosbp.polylan.h'    ! Polynomial Lanczos, only for real/Hermitian
        include 'mosbp.pist.h'       ! PIST algorithm, real/complex
        include 'mosbp.pistconv.h'
        include 'mosbp.crp.h'        ! CRP calculation, only for complex

     ! testing subroutines
        include 'mosbp.progQMR.h'
        include 'mosbp.progEig.h' 
        include 'mosbp.progLan.h' 

   end module POSBP
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
