!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Initialization subroutines for MOSB3 module            c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function initMOSB()
   integer :: ierr

   call preInit(ierr)

   if (id==rootID)   call readMOSBSTD()

   call distributeMOSB(ierr)

   call postInit(ierr)

   initMOSB = myInit() 

   if (initMOSB)   initMOSB = allocMOSB()

end function
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine finalMOSB()
   integer :: ierr

   call myFinal()

   call deallocMOSB()

   call destroyBJComm(ierr)

   call MPI_FINALIZE(ierr)

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function myInit()
   integer :: info, i
   integer :: getRecvNodesNum, getSendNodesNum 

   myInit = .false.

   allocate(sPos(nNodes,sF), ePos(nNodes,sF), bNum(nNodes,sF),   &
            locDim(nNodes,sF), stat= info)  
   if (info /= 0) return

   do i = 1, sF
      call getPosition(i,sF,sN,nNodes,sLen(i),sPos(1,i), &
                ePos(1,i),bNum(1,i),locDim(1,i))      
   end do

   myInit = .true.

end function
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine myFinal()

   if (allocated(sPos)) deallocate(sPos)   

   if (allocated(ePos)) deallocate(ePos)

   if (allocated(bNum)) deallocate(bNum)

   if (allocated(locDim)) deallocate(locDim)

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine initLayers(s1, s2)

   integer, intent(IN) :: s1, s2

   integer :: info
   integer :: getSendNodesNum, getRecvNodesNum

   sendNum = getSendNodesNum(id,nNodes,sN(s1),sLen(s1),sPos(1,s1),   &
                 ePos(1,s1),bNum(1,s1),sN(s2),sLen(s2),sPos(1,s2),   &
                 ePos(1,s2),bNum(1,s2))
   recvNum = getRecvNodesNum(id,nNodes,sN(s1),sLen(s1),sPos(1,s1),   &
                 ePos(1,s1),bNum(1,s1),sN(s2),sLen(s2),sPos(1,s2),   &
                 ePos(1,s2),bNum(1,s2))   

   allocate(nInd1(sendNum),lenInd1(sendNum),locInd1(sendNum),gInd1(sendNum), &
            nInd2(recvNum),lenInd2(recvNum),locInd2(recvNum),gInd2(recvNum), &
            req1(sendNum), req2(recvNum),reqx1(sendNum),reqx2(recvNum),      &
            gIndx1(sendNum), gIndX2(recvNum),gridInd(sendNum),stat=info)     
  
   call getSendNodesInd(id,nNodes,sN(s1),sLen(s1),sPos(1,s1),ePos(1,s1),  &
          bNum(1,s1),locDim(1,s1),sN(s2),sLen(s2),sPos(1,s2),ePos(1,s2),  &
          bNum(1,s2),locDim(1,s2),sendNum,nInd1,lenInd1,locInd1,gInd1,gridInd)   

   call getRecvNodesInd(id,nNodes,sN(s1),sLen(s1),sPos(1,s1),ePos(1,s1),  &
          bNum(1,s1),locDim(1,s1),sN(s2),sLen(s2),sPos(1,s2),ePos(1,s2),  &
          bNum(1,s2),locDim(1,s2),recvNum,nInd2,lenInd2,locInd2,gInd2) 

   gIndX1(1:sendNum) = gInd1(1:sendNum)+ XTAG_SHIFT
   gIndX2(1:recvNum) = gInd2(1:recvNum)+ XTAG_SHIFT

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine finalLayers()
   
     if (allocated(nInd1))    deallocate(nInd1)
     if (allocated(lenInd1))  deallocate(lenInd1)
     if (allocated(locInd1))  deallocate(locInd1)    
     if (allocated(gInd1))    deallocate(gInd1)
     if (allocated(req1))     deallocate(req1)
     if (allocated(reqX1))    deallocate(reqX1)
     if (allocated(gIndX1))   deallocate(gIndX1)

     if (allocated(gridInd))  deallocate(gridInd)
  
     if (allocated(nInd2))    deallocate(nInd2)
     if (allocated(lenInd2))  deallocate(lenInd2)
     if (allocated(locInd2))  deallocate(locInd2)    
     if (allocated(gInd2))    deallocate(gInd2)
     if (allocated(req2))     deallocate(req2)
     if (allocated(reqX2))    deallocate(reqX2)
     if (allocated(gIndX2))   deallocate(gIndX2)

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


