!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c          Subroutines to get the indices for Vi                  c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine mgetViColInd(nNodes,sF,sN,col,grpInd,rootInd,nodNum,blkInd,snInd)
   implicit none  
   include 'mpif.h'
   integer,intent(IN)  :: nNodes, sF, sN(sF)
   integer(kind=MPI_OFFSET_KIND), intent(IN) :: col
   integer,intent(OUT) :: grpInd(sF),rootInd(sF),nodNum(sF),blkInd(sF),snInd(sF) 

   integer(kind=MPI_OFFSET_KIND) :: gDim(sF),gBlk(sF)
   integer :: i

   gDim(1) = 1
   do i = 1, sF-1
      gDim(i+1) = gDim(i)*sN(i)
   end do

   gBlk(sF) = 1
   do i = sF, 2,-1
      gBlk(i-1) = gBlk(i)*sN(i)
   end do

   call mgetViColIndex(nNodes,sF,sN,gDim,gBlk,col,grpInd,rootInd, &
                       nodNum,blkInd,snInd)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine mgetViColIndex(nNodes, sF, sN, gDim, gBlk, col, &
                          grpInd,rootInd, nodeNum, blkInd, snInd)
   implicit none  
   include 'mpif.h'
   integer, intent(IN) :: nNodes, sF, sN(sF)
   integer(kind=MPI_OFFSET_KIND), intent(IN)  :: gDim(sF),gBlk(sF),col
   integer,intent(OUT):: grpInd(sF),rootInd(sF),nodeNum(sF),blkInd(sF),snInd(sF)
   
   integer :: i
   integer(kind=MPI_OFFSET_KIND) :: blk,pos, glen(sF) 
   integer :: pnum, qnum, cnum

   glen(1:sF) = sN(1:sF)*gDim(1:sF)
 
   do i = sF, 1, -1
      blk = (col-1)/gLen(i)+1
      pos = col - (blk-1)*gLen(i)
      snInd(i) = (pos-1)/gDim(i) + 1

      if (gBlk(i)<nNodes) then
         grpInd(i)=blk-1
         blkInd(i)=1   
       
         pnum = nNodes/gBlk(i)
         qnum = nNodes - pnum*gBlk(i)

         if (blk<=qnum) then             
             rootInd(i)=(blk-1)*(pnum+1)
             nodeNum(i)=pnum+1
         else
             rootInd(i)=qnum*(pnum+1)+(blk-qnum-1)*pnum
             nodeNum(i)=pnum
         end if
      else
         nodeNum(i)=1
         pnum = gBlk(i)/nNodes
         qnum = gBlk(i) - pnum*nNodes
         cnum = qnum * (pnum+1)

         if (blk<=cnum) then
             rootInd(i)=(blk-1)/(pnum+1)
             blkInd(i)=blk-(pnum+1)*rootInd(i)
         else
             rootInd(i)=(blk-cnum-1)/pnum+qnum
             blkInd(i)=blk-cnum-(rootInd(i)-qnum)*pnum
         end if
         grpInd(i)=rootInd(i)
      end if

   end do
end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine mgetViPosInd(nNodes,sF,sN,pos,grpInd,rootInd,nodNum,idInd,blkInd,  &
                        snInd,colInd)
   implicit none  
   include 'mpif.h'
   integer, intent(IN) :: nNodes, sF, sN(sF)
   integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
   integer, intent(OUT):: grpInd(sF),rootInd(sF),nodNum(sF),idInd(sF)
   integer, intent(OUT):: blkInd(sF),snInd(sF),colInd(sF)

   integer(kind=MPI_OFFSET_KIND) :: gDim(sF),gBlk(sF)
   integer :: i

   gDim(1) = 1
   do i = 1, sF-1
      gDim(i+1) = gDim(i)*sN(i)
   end do

   gBlk(sF) = 1
   do i = sF, 2,-1
      gBlk(i-1) = gBlk(i)*sN(i)
   end do

   call mgetViPosIndex(nNodes,sF,sN,gDim,gBlk,pos,grpInd,rootInd,nodNum,idInd, &
                       blkInd,snInd,colInd)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine mgetViPosIndex(nNodes, sF, sN, gDim, gBlk, pos, grpInd,       &
                          rootInd, nodeNum, idInd,blkInd, snInd, colInd)
   implicit none  
   include 'mpif.h'
   integer, intent(IN) :: nNodes, sF, sN(sF)
   integer(kind=MPI_OFFSET_KIND), intent(IN)  :: gDim(sF),gBlk(sF),pos
   integer, intent(OUT) :: grpInd(sF),rootInd(sF),nodeNum(sF),idInd(sF)
   integer, intent(OUT) :: blkInd(sF),snInd(sF),colInd(sF)
   
   integer :: i
   integer(kind=MPI_OFFSET_KIND) :: blk, colPos, nodPos, gLen(sF)
   integer :: pnum, qnum, cnum

   glen(1:sF) = gDim(1:sF)*sN(1:sF)

   do i = sF, 1, -1
      blk = (pos-1)/gLen(i)+1
      colPos = pos - (blk-1)*gLen(i)
      snInd(i) = (colPos-1)/gDim(i) + 1
      nodPos = colPos - (snInd(i)-1)*gDim(i)

      if (gBlk(i)<nNodes) then
         grpInd(i)=blk-1
         blkInd(i)=1   
       
         pnum = nNodes/gBlk(i)
         qnum = nNodes - pnum*gBlk(i)

         if (blk<=qnum) then
             rootInd(i)=(blk-1)*(pnum+1)
             nodeNum(i)=pnum+1
         else
             rootInd(i)=qnum*(pnum+1)+(blk-qnum-1)*pnum
             nodeNum(i)=pnum
         end if

         if (nodeNum(i)==1) then
            colInd(i)=nodPos; idInd(i)=rootInd(i)
         else
            pnum   = gDim(i)/nodeNum(i)
            qnum   = gDim(i) - pnum*nodeNum(i)
            colPos = qnum*(pnum+1)
            if (nodPos<=colPos) then
                idInd(i) = (nodpos-1)/(pnum+1)
                colInd(i) = nodPos - idInd(i)*(pnum+1)
            else
                nodPos = nodPos-colPos
                idInd(i) = (nodPos-1)/pnum+qnum
                colInd(i) = nodPos - (idInd(i)-qnum)*pnum
            end if
            idInd(i) = idInd(i)+rootInd(i) 
         end if

      else
         nodeNum(i)=1;  
         pnum = gBlk(i)/nNodes
         qnum = gBlk(i) - pnum*nNodes
         cnum = qnum * (pnum+1)

         if (blk<=cnum) then
             rootInd(i)=(blk-1)/(pnum+1)
             blkInd(i)=blk-(pnum+1)*rootInd(i)
         else
             rootInd(i)=(blk-cnum-1)/pnum+qnum
             blkInd(i)=blk-cnum-(rootInd(i)-qnum)*pnum
         end if  

         idInd(i)  = rootInd(i)
         colInd(i) = nodPos
         grpInd(i) = rootInd(i) 
      end if

   end do
end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine mgetViPos(level,nNodes, sF, sN, grpInd, rootInd, nodeNum,&
                           idInd, blkInd, snInd, colInd,pos)
   implicit none  
   include 'mpif.h'
   integer, intent(IN) :: level, nNodes, sF, sN(sF)
   integer, intent(IN) :: grpInd,rootInd,nodeNum,idInd
   integer, intent(IN) :: blkInd,snInd,colInd
   integer(kind=MPI_OFFSET_KIND), intent(OUT) :: pos

   integer(kind=MPI_OFFSET_KIND) :: gDim, gBlk
   integer :: i

   gDim = 1
   do i = 1, level-1
      gDim = gDim*sN(i)
   end do

   gBlk = 1
   do i = sF, level+1,-1
      gBlk = gBlk*sN(i)
   end do

   call mgetViPosition(nNodes, sN(level), gDim, gBlk, grpInd,rootInd, nodeNum,&
                           idInd,blkInd, snInd, colInd,pos)

end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine mgetViPosition(nNodes, sN, gDim, gBlk, grpInd, rootInd, nodeNum,   &
                           idInd,blkInd, snInd, colInd,pos)
   implicit none  
   include 'mpif.h'
   integer, intent(IN) :: nNodes, sN
   integer(kind=MPI_OFFSET_KIND), intent(IN)  :: gDim,gBlk
   integer, intent(IN) :: grpInd,rootInd,nodeNum,idInd,blkInd,snInd,colInd
   integer(kind=MPI_OFFSET_KIND), intent(OUT) :: pos

   integer(kind=MPI_OFFSET_KIND) :: glen, blk0, m0Ind   
   integer :: n0Ind
   integer :: pnum, qnum, cnum, dnum, enum
   integer(kind=MPI_OFFSET_KIND) :: pblk,qblk

   glen  = sN*gDim;       n0Ind = snInd

   if (nNodes==gBlk) then
      blk0  = rootInd+1
      m0Ind = colInd
   else
      if (gBlk<nNodes) then
          blk0 = (grpInd+1)

          pnum = nNodes/gBlk
          qnum = nNodes - pnum*gBlk
          cnum = qnum*(pnum+1)
 
          if (idInd<cnum) then
              dnum = pnum + 1
          else
              dnum = pnum
          end if
          
          if (dnum==1) then
             m0Ind = colInd
          else
             pblk = gDim/dnum
             qblk = gDim - pblk*dnum
             enum = idInd - rootInd
 
             if (enum<qblk) then
                m0Ind = enum*(pblk+1)+colInd
             else
                m0Ind = qblk*(pblk+1)+(enum-qblk)*pblk+colInd
             end if             
          end if
      else
          m0Ind = colInd
          pblk = gblk/nNodes
          qblk = gblk - pblk*nNodes
 
         if (idInd<qblk) then
            blk0 = idInd*(pblk+1)+blkInd
         else
            blk0 = qblk*(pblk+1)+(idInd-qblk)*pblk + blkInd
         end if
      end if

   end if

   pos = (blk0-1) * gLen + (n0Ind-1)*gDim + m0Ind

end 
