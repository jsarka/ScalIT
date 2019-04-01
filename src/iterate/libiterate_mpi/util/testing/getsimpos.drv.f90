!
! testing program for getsimpos
!

program test_smp
   implicit none
   include "mpif.h"
   integer, parameter  :: NMAX = 20
   integer :: sF, sN(NMAX)
   integer :: nNodes
   integer :: i, j

   integer(kind=MPI_OFFSET_KIND) :: sLen,sPos,ePos
   integer :: locDim, bNum

   read(*,*) nNodes, sF
   read(*,*) (sN(i), i=1, sF)

   print *
   print *, ' # of Layers:', sF
   print *, ' Layer configuration:'
   print *, sN(1:sF)
   print *, ' # of nodes:', nNodes
 
   do i = 1, sF
      print *
      print *, ' Layer = ', i
      print *, ' nodeID     sLen  locLen bNum     sPos     ePos '
      do j = 0, nNodes - 1
         call getSimplePos(i,sF, sN, nNodes,j,sLen,sPos,ePos,bNum, locDim)      
         write(*,100) j, sLen, locDim, bNum, sPos, ePos 
      end do

   end do

   100 format(I5,1X, I10, 1X, 2(I5,1X), 2(I10,1X))

end program
