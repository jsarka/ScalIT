!
! test program for getPos
!
program test_getpos
   include 'mpif.h'
   integer :: N, sNode 
   integer, allocatable :: sN(:)
   
   integer :: i, j, k
   integer(kind=MPI_OFFSET_KIND), allocatable :: sPos1(:), ePos1(:)
   integer(kind=MPI_OFFSET_KIND), allocatable :: sPos2(:), ePos2(:)
   integer, allocatable :: locDim1(:), locDim2(:), bNum1(:), bNum2(:)
   integer(kind=MPI_OFFSET_KIND) :: nDim1, nDim2

   integer :: n1, n2
   integer :: getRecvNodesNum, getSendNodesNum

   print *, ' Input number of degrees of freedom: N='
   read(*,*) N
   allocate(sN(N))
   print *, ' Input size for each degree of freedom: sN(1:N):'
   read(*,*) sN(1:N)
   print *, ' Input Number of nodes:'
   read(*,*) sNode
   allocate(sPos1(sNode),ePos1(sNode),bNum1(sNode),sPos2(sNode),ePos2(sNode),   &
            bNum2(sNode),locDim1(sNode),locDim2(sNode))

   do i = 1, N
       print *, '=========================================================='
       call getPosition(i, N, sN, sNode, nDim1, sPos1, ePos1, bNum1, locDim1)
       print *, ' Call getPosition:'
       j = i-1
       if (j==0) j=N
       call getPosition(j, N, sN, sNode, nDim2, sPos2, ePos2, bNum2,locDim2)
       print *, ' Configuration for layer=',i
       do k=1, sNode
          print *, nDim1, sPos1(k), ePos1(k), bNum1(k)
       end do
       print *, ' Configuration for layer=',j
       do k=1, sNode
          print *, nDim2, sPos2(k), ePos2(k), bNum2(k)
       end do

       do k = 0, sNode-1
          n1 = getSendNodesNum(k,sNode,sN(i), nDim1, sPos1, ePos1, bNum1,     &
                               sN(j),nDim2, sPos2, ePos2, bNum2)
          n2 = getRecvNodesNum(k,sNode,sN(i), nDim1, sPos1, ePos1, bNum1,     &
                               sN(j),nDim2, sPos2, ePos2, bNum2)

          print *, ' For computing node=', k
          print *, ' Communication between layer=',i, ' and layer=', j
          print *, '  Send Communication:', n1,' Receive Communication:', n2
          print * , '--------------------------------------------------------'
          print *
       end do
   end do

   print *, 'finish!'

   deallocate(sN, sPos1, ePos1, bNum1, sPos2, ePos2, bNum2, locDim1, locDim2)

end program
