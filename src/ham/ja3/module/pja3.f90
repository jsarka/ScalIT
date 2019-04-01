!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Module to deal with tri-atomic molecules          c
!c       Calculate Matrix elements: MPI version          c
!c           For MPI 1, without Parallel IO              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
module PJA3
    implicit none
    include 'mpif.h'
    private
    include 'ja3.data'   !...1
    include 'mja3.data'    

    interface readJA3
       module procedure readJA3StdIO
       module procedure readJA3File
    end interface

    public :: myid, rootID, nproc
    public :: mInit, mFinal, printJA3
    public :: mCalSaveH0, mCalSaveHGm

  contains
    include 'ja3.init'   !...2
    include 'ja3.io'     !...3
    include 'mja3.io'
    include 'mja3.h0'
    include 'ja3.h0'     !...4
    include 'ja3.h0dvr'  !...5
    include 'ja3.hgm0'   !...6
    include 'ja3.hre'    !...7
    include 'ja3.hgm'    !...8
    include 'ja3.hgmij'  !...9
    include 'ja3.spline' !...10
    include 'ja3.ap'     !...11
!    include 'pja3.hgm'
    include 'pja30.hgm'

!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
   logical function mInit()
   double precision :: db
   integer :: ierr, MPI_THREAD_MULTIPLE, iprovided
   integer :: MPI_THREAD_FUNNELED, MPI_THREAD_SERIALIZED

      call MPI_INIT(ierr)
!c    call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, iprovided, ierr)
!c    call MPI_INIT_THREAD(MPI_THREAD_FUNNELED , iprovided, ierr)
!c    call MPI_INIT_THREAD(MPI_THREAD_SERIALIZED , iprovided, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid,  ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

      inquire(IOLENGTH=dbSize) db
!c    call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, dbSize, ierr)

      select case (nproc)
      case (1)
         workID(1:2)=rootID
      case (2)
         workID(1)=rootID; workID(2)=rootID+1
      case default
         workID(1)=rootID; workID(2)=rootID+1
      end select

      if (rootID==myid)   call readJA3()
      call DistributeInput()

      mInit = allocJA3()

      if ((useSP).AND.(mInit)) then
         if (myid==rootID) then
             mInit = readSpVlr()
             if (mInit) mInit=readSpVBr()
         end if
      
         call MPI_BCAST(mInit,1,MPI_LOGICAL,rootID,MPI_COMM_WORLD,ierr)
      
         if (.NOT. mInit) return
         call DistributeDVRData()
      end if
   end function
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
   subroutine mFinal()

      call deallocJA3()
      if (allocated(REVmat))  deallocate(REVmat)
      call MPI_FINALIZE(ierr)

   end subroutine
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
end module
