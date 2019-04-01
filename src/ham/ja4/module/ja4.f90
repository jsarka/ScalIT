!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Module to deal with tetra-atomic molecules        c
!c            Calculate Matrix elements                  c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
module JA4
  implicit none
  private

  include 'ja4.data'    ! data for JA4 module

  interface readJA4       ! Read input parameters
       module procedure readJA4StdIO
       module procedure readJA4File
  end interface

  public :: initJA4, finalJA4, printJA4
  public :: calSaveH0, calSaveHGm
  
  contains

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    include 'ja4.init'     ! allocate memory, and initialize data
    include 'ja4.io'       ! read input parameters and spline data

    include 'ja4.h0'       ! interface for H0 calculation
    include 'ja4.h0dvr'    ! perform DVR calculation
    include 'ja4.spline'   ! spline functions

    include 'ja4.hgm0'     ! Interface for HGM calculation
    include 'ja4.hgm'      ! calculate HGM at each point
    include 'ja4.hre'
    include 'ja4.hgmij'    ! do the real work; use more momory, but faster

 !   include 'ja4.hgm1'    ! use less memory but slower 
 !   include 'ja4.hgm4'    ! same as ja4.hgm2.h but calculate directly 
                             ! used only for testing
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    logical function initJA4()

        call readJA4()

        initJA4 = allocJA4()

        if (useSP)  then
           if ((.NOT.fixR1).AND. initJA4) initJA4 = readSpVlr1()
           if ((.NOT.fixR2).AND. initJA4) initJA4 = readSpVlr2()
           if ((.NOT.fixBR).AND. initJA4) initJA4 = readSpVBr()
        end if

    end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine finalJA4()

        call deallocJA4()

    end subroutine

end module
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
