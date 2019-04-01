!cccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Define some commonly used data structures    c
!cccccccccccccccccccccccccccccccccccccccccccccccccc

module OSBType
   implicit none
   integer, parameter :: FMAX = 20
   integer, parameter :: MAX_FNAME = 128
   integer, parameter :: STDFH =5

   TYPE CMemInfo
       integer :: mStart, mEnd, mSize
   END TYPE 

   TYPE CConvSimple 
      integer :: mMax
      double precision :: mTol
   END TYPE

   TYPE CConv
       double precision :: mE0, mTol
       integer :: mStart, mStep, mMax
       integer :: mNum, mGap
   END TYPE

   TYPE COsbw
      double precision :: mE0, mDE, mBeta
      integer :: mCnt
   END TYPE
 
   TYPE CDataInfo
      integer :: mLevel
      integer :: mSize(FMAX), mStart(FMAX), mEnd(FMAX)
      integer :: mLen, mMaxSize
   END TYPE

  contains
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calDataInfo(level, dInfo)
      integer, intent(IN) :: level
      TYPE(CDataInfo), intent(INOUT) :: dInfo

      integer :: i

      dInfo%mLevel    = level
      dInfo%mStart(1) = 1
      do i=1, level-1
         dInfo%mStart(i+1) = dInfo%mStart(i) + dInfo%mSize(i)
      end do

      dInfo%mEnd(1:level) = dInfo%mStart(1:level) + dInfo%mSize(1:level)-1
      dInfo%mLen     = dInfo%mEnd(level)
      dInfo%mMaxSize = MaxVal(dInfo%mSize(1:level))

    end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine printDataInfo(dInfo)
       TYPE(CDataInfo), intent(IN) :: dInfo

       integer :: i

       print *, ' ============================================== '
       print *, '            Layer Information of Data      '
       print *, ' ==============================================='
       print *, "    Number of Data Layers: ", dInfo%mLevel
       print *, "    Total Length of Data:", dInfo%mLen
       print *, "    Max. Length for 1 Layer:", dInfo%mMaxSize
       print *, " -----------------------------------------------"
       print *, "      Layer      Start       End       Size "
       print *, " -----------------------------------------------"
       do i = 1, dInfo%mLevel
          write(*,100) i, dInfo%mStart(i), dInfo%mEnd(i), dInfo%mSize(i)
       end do
       print *, ' ==========  Finish Layer Info  ================'

      100 FORMAT(2X, I6, 2X, 3(I10,1X))

    end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine printConv(conv)
       TYPE(CConv), intent(IN) :: conv

       print *, '=================================== '
       print *, '     Convergence Parameters         '
       print *, '===================================='
       print *
       write(*,100) conv%mStart, conv%mStep, conv%mMax
       write(*,110) conv%mNum, conv%mGap, conv%mTol
       print *

       100 FORMAT( " Start #: ", I6, " Step #:", I6, &
                    " Max.  #:" ,I6)
       110 FORMAT( " # of Compared Values:", I6,   &
                   " Gap between two sets to be compared:", &
                   "  Conv. Tol: ", F10.6)

    end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine printSimConv(conv)
       TYPE(CConvSimple), intent(IN) :: conv

       print *, '=================================== '
       print *, '   Simple Convergence Parameters    '
       print *, '===================================='
       print *

       write(*,100) conv%mMax, conv%mTol
       print *

       100 FORMAT( " Max. # of iteration: ", I10, " Conv. Tol.:", F10.6)
    end subroutine

end module 

