!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c          Get the Positions of the vector at each layer              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getPosition(level, sF, sN, nNode, sLen, sPos, ePos, bNum, locDim)
    implicit none
    include 'mpif.h'
    integer, intent(IN) :: level, sF, sN(sF), nNode
    integer(kind=MPI_OFFSET_KIND),intent(OUT) :: sLen,sPos(nNode),ePos(nNode)
    integer, intent(OUT) :: bNum(nNode), locDim(nNode)    

    integer :: i, j
    integer :: pNum, qNum, cNum, nSize, mSize
    integer(kind=MPI_OFFSET_KIND) :: nBlk

    sLen = 1;  nBlk = 1
    do i = 1, level-1
       sLen = sLen*sN(i)
    end do
    
    do i = sF, level+1, -1
       nBlk = nBlk*sN(i)
    end do

    if (nBlk>=nNode) then   ! each node has more than one node
        pNum = nBlk/nNode;     qNum = nBlk-pNum*nNode
        bNum(1:nNode)=pNum;  locDim(1:nNode) = sLen
  
        do i = 1, qNum
           bNum(i) = pNum + 1
        end do
 
        sPos(1)=1; ePos(1)=bNum(1)*sLen*sN(level)
        do i = 2, nNode
           sPos(i)=ePos(i-1)+1
           ePos(i) = sPos(i)+bNum(i)*sLen*sN(level)-1
        end do

    else     ! more than one nodes deal with one block
        pNum = nNode/nBlk;  qNum = nNode-pNum*nBlk

        bNum(1:nNode)=0; sPos(1)=1; 

        cNum = 1;  pNum = pNum + 1
        nSize = sLen/pNum; mSize = sLen-nSize*pNum

        do i = 1, qNum   ! for each block
           do j = 1, mSize
              ePos(cNum)=sPos(cNum)+nSize
              locDim(cNum)=nSize+1
              cNum=cNum+1
              if (cNum>nNode) return
              sPos(cNum)=ePos(cNum-1)+1
           end do

           do j = mSize+1, pNum
              ePos(cNum)=sPos(cNum)+nSize-1
              locDim(cNum)=nSize
              cNum = cNum + 1
              if (cNum>nNode) return
              sPos(cNum)=ePos(cNum-1)+1
           end do
           sPos(cNum) = ePos(cNum-1) + (sN(level)-1)*sLen+1
        end do
        
        pNum = nNode/nBlk
        nSize = sLen/pNum; mSize = sLen-nSize*pNum

        do i = qNum+1, nBlk   ! for each block
           if (pNum==1) then
              ePos(cNum)=sPos(cNum)+sLen*sN(level)-1
              locDim(cNum)=sLen;     bNum(cNum)=1
              cNum = cNum + 1
              if (cNum>nNode) return
              sPos(cNum)=ePos(cNum-1)+1
           else
              do j = 1, mSize
                 ePos(cNum)=sPos(cNum)+nSize
                 locDim(cNum)=nSize+1
                 cNum = cNum + 1
                 if (cNum>nNode) return
                 sPos(cNum)=ePos(cNum-1)+1
              end do

              do j = mSize+1, pNum
                 ePos(cNum)=sPos(cNum)+nSize-1
                 locDim(cNum)=nSize
                 cNum = cNum + 1
                 if (cNum>nNode) return
                 sPos(cNum)=ePos(cNum-1)+1
              end do
              sPos(cNum) = ePos(cNum-1) + (sN(level)-1)*sLen+1
           end if

        end do
    end if
end 

