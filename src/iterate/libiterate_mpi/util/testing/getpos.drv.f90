!
! test program for getPos
!
program test_getpos
   include 'mpif.h'
   integer :: N, sNode
   
   integer :: i, bnum1, mDim1
   integer(kind=MPI_OFFSET_KIND), allocatable :: sPos(:), ePos(:)
   integer, allocatable :: sN(:), mDim(:), bNum(:)
   integer(kind=MPI_OFFSET_KIND) :: nDim, sPos1, ePos1, nDim1


 !  print *, ' Input number of the degrees of freedom:'
   read(*,*) sNode, N
   allocate(sN(N)) 
   read(*,*) sN(1:N)

   print *
   print *, ' # of Layers:', N
   print *, ' Layer configuration:'
   print *, sN(1:N)
   print *, ' # of nodes:', sNode

   allocate(sPos(sNode), ePos(sNode), bNum(sNode), mDim(sNode))

   do i = 1, N
       call getPosition(i, N, sN, sNode, nDim, sPos, ePos, bNum, mDim)
       print *
       print *, ' Configuration for layer=',i
!       print *, ' gDim   start ,   end ,   block Number  locDim '
       do j=1, sNode
          call getSimplePos(i,N,sN,sNode,j-1,nDim1,sPos1,ePos1,bNum1,mDim1)
          if (((nDim1-nDim)/=0) .or. ((sPos1-sPos(j))/=0) .or.     &
               ((ePos1-ePos(j))/=0) .or. ((bNum(j)-bnum1)/=0) .or. &
               ((mDim1-mDim(j))/=0) ) then
               print *, ' Error in getPosition subroutine'
          else
               print *, ' OK:layer=',i,' id=', j-1
          end if
          !print *, nDim, sPos(j), ePos(j), bNum(j), mDim(j)

       end do
       
   end do

   print *, 'finish!'

   deallocate(sPos, ePos, bNum, mDim, sN)

end program
