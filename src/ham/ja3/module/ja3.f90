!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c              Module to deal with tri-atomic molecules             c
!c                     Calculate Matrix elements                     c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

module JA3
    implicit none
    private 
    include 'ja3.data'   !...1  Definition of variables 

    interface readJA3
       module procedure readJA3StdIO
       module procedure readJA3File
    end interface

    public :: initJA3, finalJA3, printJA3 
    public :: calSaveH0,  calSaveHGM

  contains 
    include 'ja3.init'     !...2 allocate/deallocate memory
    include 'ja3.io'       !...3 read/write data
    include 'ja3.h0'       !...4 saveHlr, saveHBR interface
    include 'ja3.h0dvr'    !...5 do DVR computation
    include 'ja3.hgm0'     !...6 saveHGm interface
    include 'ja3.hre'      !...7
    include 'ja3.hgm'      !...8 interfaces for Hgm calculation
    include 'ja3.hgmij'    !...9 do the real work
    include 'ja3.spline'   !...10 related to spline function
    include 'ja3.ap'       !...11 calculate absorption potential
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    logical function initJA3()

       call readJA3()
       initJA3=allocJA3()
       if ((useSP) .AND. (initJA3) )  then
          initJA3 = readSpVlr()
          if (initJA3) initJA3=readSpVBr()
       end if

    end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine finalJA3()

       call deallocJA3()

    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

end module
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
