!
! Testing mindex.f90
!
program test_mindex
   use mosbtype
   implicit none
   include 'mpif.h'

   integer :: nNodes, sF, sN(FMAX)
   integer  :: i, j
   integer(kind=MPI_OFFSET_KIND) :: pos, pos1
   TYPE(GDataInfo) :: myData
   integer :: ind(FMAX)

   read(*,*) nNodes, sF
   read(*,*) sN(1:sF)

   call calGData(sF,sN,myData)
   call printGData(myData)

   do j = 1, sF
      print *, ' Layer=',j
      print *, '-----------------------'
      print *, '  Pos       Ind        '
      print *, '-----------------------'
      do i = 1, myData%gN
         pos = i
         call mgetVxColInd(sF,sN,pos,ind)     
         write(*,10) pos,ind(j)
      end do
   end do

   10 format (I15, I10)

end
