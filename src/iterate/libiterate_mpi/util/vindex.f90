!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c          Subroutines to get the indices for Vi                  c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getViColInd(sF, sN, col, blkInd, snInd)
   implicit none  
   include 'mpif.h'
   integer, intent(IN) :: sF, sN(sF)
   integer(kind=MPI_OFFSET_KIND), intent(IN)  :: col
   integer(kind=MPI_OFFSET_KIND), intent(OUT) :: blkInd(1)
   integer, intent(OUT):: snInd(sF)

   integer(kind=MPI_OFFSET_KIND) :: gDim(sF)
   integer :: i

   gDim(1)=1
   do i = 1, sF-1
      gDim(i+1)=gDim(i)*sN(i)
   end do

   call getViColIndex(sF,sN,gDim,col,blkInd,snInd)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getViColIndex(sF, sN, gDim, col, blkInd,  snInd)
   implicit none  
   include 'mpif.h'
   integer, intent(IN) :: sF, sN(sF)
   integer(kind=MPI_OFFSET_KIND), intent(IN)  :: gDim(sF), col
   integer(kind=MPI_OFFSET_KIND), intent(OUT) :: blkInd(sF)
   integer, intent(OUT):: snInd(sF)
   
   integer :: i
   integer(kind=MPI_OFFSET_KIND) :: pos, glen(sF)

   glen(1:sF) = gDim(1:sF)*sN(1:sF)

   do i = sF, 1, -1
      blkInd(i) = (col-1)/glen(i)+1
      pos = col - (blkInd(i)-1)*glen(i) 
      snInd(i) =  ((pos-1)/gDim(i)) + 1    
   end do
end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getViPosInd(sF, sN, col, blkInd, sNInd, colInd)
   implicit none  
   include 'mpif.h'
   integer, intent(IN) :: sF, sN(sF)
   integer(kind=MPI_OFFSET_KIND), intent(IN) :: col
   integer, intent(OUT) :: sNInd(sF)
   integer(kind=MPI_OFFSET_KIND), intent(OUT):: blkInd(sF), colInd(sF)

   integer(kind=MPI_OFFSET_KIND) :: gDim(sF)
   integer :: i

   gDim(1)=1
   do i = 1, sF-1
      gDim(i+1)=gDim(i)*sN(i)
   end do

   call getViPosIndex(sF,sN,gDim,col,blkInd,sNInd,colInd)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getViPosIndex(sF, sN, gDim, col, blkInd, sNInd, colInd)
   implicit none  
   include 'mpif.h'
   integer, intent(IN) :: sF, sN(sF)
   integer(kind=MPI_OFFSET_KIND), intent(IN)  :: gDim(sF), col
   integer(kind=MPI_OFFSET_KIND), intent(OUT) :: colInd(sF),blkInd(sF)
   integer, intent(OUT):: snInd(sF)
   
   integer :: i
   integer(kind=MPI_OFFSET_KIND) :: pos,glen(sF)

   glen(1:sF) = gDim(1:sF)*sN(1:sF)

   do i = sF, 1, -1
      blkInd(i) = (col-1)/glen(i)+1
      pos = col - (blkInd(i)-1)*glen(i)  
      snInd(i) = (pos-1)/gDim(i) + 1
      colInd(i) = pos - (snInd(i)-1)*gDim(i)   
   end do
end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getViPos(sF, sN, level, blkInd, snInd, colInd, pos)
   implicit none  
   include 'mpif.h'
   integer, intent(IN) :: sF, sN(sF), level, snInd
   integer(kind=MPI_OFFSET_KIND), intent(IN)  :: blkInd,colInd
   integer(kind=MPI_OFFSET_KIND), intent(OUT) :: pos
   
   integer :: i

   integer(kind=MPI_OFFSET_KIND) :: gDim

   gDim = 1
   do i = 1, level-1
     gDim = gDim*sN(i)
   end do

   pos = (blkInd-1)*sN(level)*gDim + (snInd-1)*gDim + colInd

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getViPosition(sN, gDim, blkInd, snInd, colInd, pos)
   implicit none  
   include 'mpif.h'
   integer, intent(IN) :: sN, sNInd
   integer(kind=MPI_OFFSET_KIND), intent(IN)  :: gDim,colInd,blkInd
   integer(kind=MPI_OFFSET_KIND), intent(OUT) :: pos

   pos = (blkInd-1)*sN*gDim + (snInd-1)*gDim + colInd
end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
