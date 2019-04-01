!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   The grouping of nodes for each block           c
!c     This only works when node > blk              c
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine indexNode(blk, node, myid, nodeId, grpId, grpSize)
   implicit none
   integer, intent(IN)  :: blk, node, myid
   integer, intent(OUT) :: nodeID, grpID, grpSize

   integer :: num1, num2, num3, num4

   num1 = node/blk;  num2 = node - blk * num1

   if (num2==0) then
      grpID = myID/num1;  nodeID = myID-grpID*num1
      grpSize = num1
   else
      num3=num1+1; num4=num2-1
      grpID=myID/num3

      if (grpID<num2) then
         nodeID=myID-num3*grpID
         grpSize = num3
      else
         grpSize = num1
         num4   = myID - num3*num2
         nodeID = num4 - num4/num1*num1
         grpID  = num4/num1+num2       
      end if
   end if
end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine indexNodeLen(gLen, node, myid, pLen, gStart)
   implicit none
   include 'mpif.h'
   integer(kind=MPI_OFFSET_KIND),intent(IN) :: gLen
   integer, intent(IN)  :: node, myID
   integer, intent(OUT) :: pLen
   integer(kind=MPI_OFFSET_KIND),intent(OUT)::gStart

   integer :: num1,num2,num3,num4, grpID

   num1 = gLen/node;  num2 = gLen - node * num1

   if (num2==0) then
      pLen = num1;  gStart = myID*pLen + 1
   else
      if (myID<num2) then
         pLen=num1+1
         gStart=myID*pLen+1
      else
         pLen = num1
         gStart = (num1+1)*num2+num1*(myID-num2)+1
      end if
   end if
end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getMPI1DIndex(total, nproc, myid, sStart, sEnd, sLen)
   implicit none
   integer,intent(IN)  :: total, nproc, myID
   integer,intent(OUT) :: sStart, sEnd, sLen

   integer :: num1,num2,num3,num4, grpID

   if (myid>=nproc)  return

   if (nproc==1) then
      sStart=1; sLen=total; sEnd=total; return
   end if

   num1 = total/nproc;  num2 = total - nproc * num1

   if (num2==0) then
      sLen = num1;  sStart = myID*sLen + 1
   else
      if (myID<num2) then
         sLen=num1+1;  sStart=myID*sLen+1
      else
         sLen = num1
         sStart = (num1+1)*num2+num1*(myID-num2)+1
      end if
   end if
   sEnd = sStart+sLen-1
end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getMPI1DIndexLong(gLen, nproc, myid, gStart, gEnd, pLen )
   implicit none
   include 'mpif.h'
   integer(kind=MPI_OFFSET_KIND),intent(IN) :: gLen
   integer, intent(IN)  :: nproc, myID
   integer, intent(OUT) :: pLen
   integer(kind=MPI_OFFSET_KIND),intent(OUT)::gStart, gEnd

   integer :: num1,num2,num3,num4, grpID

   num1 = gLen/nproc;  num2 = gLen - nproc * num1

   if (num2==0) then
      pLen = num1;  gStart = myID*pLen + 1
   else
      if (myID<num2) then
         pLen=num1+1
         gStart=myID*pLen+1
      else
         pLen = num1
         gStart = (num1+1)*num2+num1*(myID-num2)+1
      end if
   end if
   gEnd=gStart+pLen-1
end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

