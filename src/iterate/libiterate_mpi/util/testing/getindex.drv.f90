!cccccccccccccccccccccccccccccccccccccc
!c   test program for get***Index     c
!cccccccccccccccccccccccccccccccccccccc
program test_getindex
   implicit none
   include 'mpif.h'
   integer :: N, sNode 

   integer, allocatable :: sN(:)      
   integer(kind=MPI_OFFSET_KIND),allocatable :: sPos1(:),ePos1(:),bNum1(:)
   integer(kind=MPI_OFFSET_KIND),allocatable :: sPos2(:),ePos2(:),bNum2(:)
   integer(kind=MPI_OFFSET_KIND) :: nLen1, nLen2
   integer, allocatable :: locDim1(:), locDim2(:)
   integer, allocatable :: nnInd(:), lenInd(:), locInd(:), gridInd(:)

   integer :: n1, n2, i, j, k, kk
   integer :: getRecvNodesNum, getSendNodesNum

   print *, ' Input number of degrees of freedom: N='
   read(*,*) N
   allocate(sN(N))
   print *, ' Input size for each degree of freedom: sN(1:N):'
   read(*,*) sN(1:N)
   print *, ' Input Number of nodes:'
   read(*,*) sNode

   allocate(locDim1(sNode), locDim2(sNode),sPos1(sNode),ePos1(sNode),     &
            bNum1(sNode),sPos2(sNode),ePos2(sNode),bNum2(sNode))

!   do i = 1, N
       i = 3;
       print *, '=========================================================='
       call getPosition(i, N, sN, sNode, nLen1, sPos1, ePos1, bNum1, locDim1)
       j = i-1; j = 1
       if (j==0) j=N
       call getPosition(j, N, sN, sNode, nLen2, sPos2, ePos2, bNum2, locDim2)
       print *, ' Configuration for layer=',i
       print *, ' sLen    sPos   ePos   bNum     locDim'
       do k=1, sNode
          print *, nLen1, sPos1(k), ePos1(k), bNum1(k), locDim1(k) 
       end do
       print *, ' Configuration for layer=',j
       print *, ' sLen    sPos   ePos   bNum     locDim'
       do k=1, sNode
          print *, nLen2, sPos2(k), ePos2(k), bNum2(k), locDim2(k)
       end do

       do k = 0, sNode-1
          n1 = getSendNodesNum(k,sNode,sN(i), nLen1, sPos1, ePos1, bNum1,  &
                               sN(j),nLen2, sPos2, ePos2, bNum2)
          allocate(nnInd(n1),lenInd(n1),locInd(n1),gridInd(n1))
          call getSendNodesInd(k,sNode,sN(i), nLen1, sPos1, ePos1, bNum1, locDim1,   &
                               sN(j),nLen2, sPos2, ePos2, bNum2, locDim2, n1, nnInd, &
                               lenInd, locInd, gridInd)

          print *, ' For computing node=', k
          print *, ' Initial data layer = ',i, '. Final data layer = ', j
          print *, ' Nodes communication involved: N1=', n1
          print *, ' Initila layer = ', i, ' Final Layer = ', j
          print *, ' Computing node = ', k , ' send data to the following nodes.'
          print *, ' Targ. Node    Length   StartPos  GridTag'
          do kk=1, n1
             print *, nnInd(kk), lenInd(kk), locInd(kk), gridInd(kk)
          end do
          deallocate(nnInd, lenInd, locInd, gridInd)

          n2 = getRecvNodesNum(k,sNode,sN(i), nLen1, sPos1, ePos1, bNum1,      &
                               sN(j),nLen2, sPos2, ePos2, bNum2)
          allocate(nnInd(n2),lenInd(n2),locInd(n2),gridInd(n2))
          call getRecvNodesInd(k,sNode,sN(i), nLen1, sPos1, ePos1, bNum1,locDim1,    &
                               sN(j),nLen2, sPos2, ePos2, bNum2, locDim2, n2, nnInd, &
                               lenInd, locInd, gridInd)

          print *
          print *, ' For computing node=', k
          print *, ' Initial data layer = ',i, '. Final data layer = ', j
          print *, ' Nodes communication involved: N2=', n2
          print *, ' Initila layer = ', i, ' Final Layer = ', j
          print *, ' Computing node = ', k , ' receive data from the following nodes.'
          print *, ' Targ. Node    Length   StartPos  GridTag'
          do kk=1, n2
             print *, nnInd(kk), lenInd(kk), locInd(kk), gridInd(kk)
          end do

          deallocate(nnInd, lenInd, locInd, gridInd)

          print * , '--------------------------------------------------------'
          print *
       end do
!   end do

   print *, 'finish!'

   deallocate(locDim1, locDim2, sN, sPos1, ePos1, bNum1, sPos2, ePos2, bNum2)


end program
