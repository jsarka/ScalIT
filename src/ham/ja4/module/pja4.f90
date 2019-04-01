!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Module to deal with tetra-atomic molecules        c
!c       Calculate Matrix elements:MPI version           c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
module PJA4
    implicit none
    include 'mpif.h'
    private
    include 'ja4.data'
    include 'mja4.data'    ! data for JA4 module

    interface readMJA4
       module procedure readJA4StdIO
       module procedure readJA4File
    end interface

    public :: myid, rootID, nproc
    public :: mInit, mFinal, printJA4
    public :: MCalSaveH0, MCalSaveHGm

  contains
    include 'ja4.init'
    include 'ja4.io'
    include 'mja4.io'      ! read/write data   
    include 'ja4.h0dvr'
    include 'ja4.spline'
    include 'ja4.h0'
    include 'mja4.hre'
    include 'ja4.hgm'
    include 'ja4.hgmij'
    include 'ja4.hgm0'
    include 'mja4.h0'
!    include 'pja4.hgm'
    include 'pja40.hgm'

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  logical function mInit()
       double precision :: db

       call MPI_Init(ierr)
       call MPI_Comm_Rank(MPI_COMM_WORLD, myid,  ierr)
       call MPI_Comm_Size(MPI_COMM_WORLD, nproc, ierr)

       inquire(IOLENGTH=dbSize) db
!       call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, dbSize, ierr)
       
       select case (nproc)
       case (1)
          workID(1:2)=0
       case (2)
          workID(1)=0; workID(2)=1     
       case default
          workID(1)=0; workID(2)=1
       end select

       if (myid==rootID)  call readMJA4() 
       call DistributeInput()
       
       mInit = allocJA4()

       if (useSP .AND. mInit) then
          if (myid==rootID) then
             mInit = readSpVlr1()
             if (mInit) mInit=readSpVlr2()
             if (mInit) mInit=readSpVBr()
          end if
          call MPI_BCAST(mInit,1,MPI_LOGICAL,rootID,MPI_COMM_WORLD,ierr)
          if (.NOT. mInit) return
          call DistributeDVRData()
       end if

    end function
!*************************************************************************

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine mFinal()

       call deallocJA4()

       if (allocated(REVmat))  deallocate(REVmat)

       call MPI_Finalize(ierr)

    end subroutine

end module
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


