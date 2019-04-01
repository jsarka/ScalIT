!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    OSB implementation using direct send/recv      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
! below change from 2**28 limits to 4096 cores -CP Mar 2016
  integer, parameter :: XTAG_SHIFT = 2**12

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccc  Number of nodes involved in send   cccccccccc 
!cccccc   and receive data at each layer    cccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc

  integer(KIND=MPI_OFFSET_KIND) :: sLen(FMAX)

  integer(KIND=MPI_OFFSET_KIND),allocatable :: sPos(:,:),ePos(:,:)

  integer, allocatable :: bNum(:,:), locDim(:,:)  
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   integer,allocatable :: nInd1(:),lenInd1(:),locInd1(:),gInd1(:),   &
                          req1(:),reqx1(:),gIndx1(:), gridInd(:)
   integer,allocatable :: nInd2(:),lenInd2(:),locInd2(:),gInd2(:),   &
                          req2(:),reqx2(:),gIndx2(:)

   integer :: sendNum, recvNum
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

