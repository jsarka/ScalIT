!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                Get indices for V*X operations                        c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine mgetVxColInd(sF,sN,col,ind)
   implicit none  
   include 'mpif.h'
   integer,intent(IN)  ::  sF, sN(sF)
   integer(kind=MPI_OFFSET_KIND), intent(IN) :: col
   integer,intent(OUT) :: Ind(sF) 

   integer(kind=MPI_OFFSET_KIND) :: gDim(sF)
   integer :: i

   gDim(1) = 1
   do i = 1, sF-1
      gDim(i+1) = gDim(i)*sN(i)
   end do  

   call mgetVxColIndex(sF,sN,gDim,col,ind)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine mgetVxColIndex(sF, sN, gDim, col, ind)
   implicit none  
   include 'mpif.h'
   integer, intent(IN) :: sF, sN(sF)
   integer(kind=MPI_OFFSET_KIND), intent(IN)  :: gDim(sF),col
   integer,intent(OUT):: ind(sF)
   
   integer(kind=MPI_OFFSET_KIND)  :: gLen(sF)

   gLen(1:sF) = gDim(1:sF)*sN(1:sF)

   call mgetVxColIndices(sF,sN,gDim,gLen,col, ind)

end 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine mgetVxColIndices(sF, sN, gDim, gLen, col, ind)
   implicit none  
   include 'mpif.h'
   integer, intent(IN) :: sF, sN(sF)
   integer(kind=MPI_OFFSET_KIND), intent(IN)  :: gDim(sF),gLen(sF),col
   integer,intent(OUT):: ind(sF)

   integer :: i
   integer(kind=MPI_OFFSET_KIND) :: blk, pos

   do i = 1, sF
      blk = (col-1)/gLen(i)
      pos = col - blk*gLen(i)
      ind(i) = (pos-1)/gDim(i)+1
   end do

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

