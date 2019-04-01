!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c               Parameters and data for OSBW                        c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   integer,parameter :: B0Tag = 2**8    ! max Layer=2**8-1 for getVi

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ! global and local size of window states
   integer :: hwLen, phwLen    

       ! matrix elements of window states, real version
   double precision, allocatable :: HW(:,:), HWR(:,:), HWTMP(:,:)
       ! vector elements for window states, real version
   double precision, allocatable :: WX(:), WY(:)

       ! matrix elements of window states, complex version
   double complex,  allocatable :: HWCX(:,:),HWRCX(:,:),HWTMPCX(:,:)
       ! vector elements for window states, complex version
   double complex,  allocatable ::WXCX(:),WYCX(:)

       ! indices for decomposition used to calculate HW^-1
   integer, allocatable  :: WIPIV(:), WIND(:)

       ! number of windows states, and their summary order for each node
   integer, allocatable :: gCnt(:), sCnt(:), nodInd(:)

        ! local indices for all and local window states
   integer, allocatable :: pInd(:), gInd(:)  
        ! global indices for all and window states
   integer(kind=MPI_OFFSET_KIND),allocatable :: gKindex(:)
        ! local indices for local window states
   integer,allocatable :: gPIndex(:)

        ! global APP, APR and AP data
   double precision, allocatable :: gAPP(:), gAPR(:), gAP(:), gEig0(:)  

        ! indices to get Vi vector
   integer, allocatable :: grpInd(:,:),rootInd(:,:),nodeNum(:,:)
   integer, allocatable :: blkInd(:,:),sNInd(:,:)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
