!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                   Data Structures used for MOSB                        c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   TYPE :: CConvSimple
       integer :: mMax
       double precision :: mTol
   end TYPE
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   TYPE :: CConv
       double precision :: mE0, mTol
       integer :: mStart, mStep, mMax
       integer :: mNum, mGap
   end TYPE
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   TYPE :: COsbw
       double precision :: mE0, mDe, mBeta
       integer :: mCnt
   end TYPE
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   TYPE :: CSeqInfo
       integer :: mLevel
       integer :: mSize(FMAX), mStart(FMAX), mEnd(FMAX)
       integer :: mLen, mMaxSize
   end Type
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   TYPE :: GDataInfo
        integer :: sF, sN(FMAX)
        integer(kind=MPI_OFFSET_KIND) :: gDim(FMAX), gBlk(FMAX)
        integer(kind=MPI_OFFSET_KIND) :: gLen(FMAX), gN, gMaxLen
   end Type
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   TYPE :: MNodeInfo
        integer :: id, nNodes, spLevel  
        logical :: lbFlag          
        integer :: myID(FMAX), grpID(FMAX), commID(FMAX)   
        integer :: nGroup(FMAX), nodNum(FMAX), nodIDStart(FMAX)
   end TYPE
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   TYPE :: MDataInfo
        integer :: pBlk(FMAX), pDim(FMAX), pLen(FMAX), pMaxLen
        integer(kind=MPI_OFFSET_KIND) :: gBStart(FMAX), gBEnd(FMAX)    
        integer(kind=MPI_OFFSET_KIND) :: gDStart(FMAX), gDEnd(FMAX)  
   end Type
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   TYPE :: MGridInfo
        integer :: pSize(FMAX), pStart(FMAX), pEnd(FMAX), pLen, pMaxSize
        integer(kind=MPI_OFFSET_KIND) :: gStart(FMAX),gEnd(FMAX),gLen
        integer(kind=MPI_OFFSET_KIND) :: gSize(FMAX),gPos(FMAX), gMaxSize
   end Type
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

