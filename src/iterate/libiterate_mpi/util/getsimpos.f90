!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c          Get the Positions of the vector at each layer              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getSimplePos(level, sF, sN, nNode, id, sLen, sPos, ePos, bNum, locDim)
    implicit none
    include 'mpif.h'
    integer, intent(IN) :: level, sF, sN(sF), nNode, id
    integer(kind=MPI_OFFSET_KIND),intent(OUT) :: sLen,sPos,ePos
    integer, intent(OUT) :: locDim, bNum

    integer :: i, j
    integer :: pNum, qNum, cNum, dNum
    integer(kind=MPI_OFFSET_KIND) :: nBlk, pBlk

    sLen = 1;  nBlk = 1
    do i = 1, level-1
       sLen = sLen*sN(i)
    end do
    
    do i = sF, level+1, -1
       nBlk = nBlk*sN(i)
    end do

    if (nBlk>=nNode) then   ! each node has more than one block
        pNum = nBlk/nNode
        qNum = nBlk-pNum*nNode
       
        if (id < qnum) then
            bNum = pNum + 1
            pBlk = id * (pNum+1) + 1   
        else
            bNum = pNum 
            pBlk = qNum * (pnum+1) + (id-qnum)*pnum + 1
        end if 

        locDim = sLen
        sPos = (pBlk-1)*sLen*sN(level) + 1
        ePos = sPos + bNum * sLen*sN(level) -1  

    else     ! more than one nodes deal with one block

        pNum = nNode/nBlk;  qNum = nNode - pNum * nBlk              

        if ( id < qNum*(pNum+1) ) then
           pNum = pNum + 1
           cNum = id / pnum 
           dNum = id - cNum * pnum
        else
           cNum = (id-qnum*(pNum+1))/pnum + qnum
           dNum = id - qnum*(pNum+1) - (cNum-qnum)*pNum 
        end if

        if (pNum==1) then
           bNum = 1; sPos = 1
           locDim = sLen    
           sPos = cNum*sLen*sN(level) + sPos
           ePos = sPos + locDim*sN(level) - 1     
        else
           bNum = 0
           call seqDataLenLong(sLen, pnum, dnum, locDim, sPos)
           sPos = cNum*sLen*sN(level) + sPos
           ePos = sPos + locDim - 1 
        end if        

    end if
end 

