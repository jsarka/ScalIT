!ccccccccccccccccccccccccccccccccccccccccccccccc
!c        Testing program for mosbtype         c
!ccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccc
!c   comment out MPI_COMM_SPLIT function call   c
!c    at line 49 in file mosbtype.cal.h         c
!cccccccccccccccccccccccccccccccccccccccccccccccc
program mosb_seq
   use mosbtype
   implicit none

   integer :: opt, i, j

   logical :: sSt
   integer :: nNodes, sF, sN(FMAX)
   Type(GDataInfo) :: myGd
   TYPE(MNodeInfo) :: myNode
   TYPE(MDataInfo) :: myData
   TYPE(MGridInfo) :: myHosb

   print *,'  =================================='
   print *,'        Testing MOSBType:calVOSB   '
   print *, ' =================================='
   print *
   print *, ' Read input parameters from standard input. '

   read(*,*) nNodes, sF
   read(*,*) sN(1:sF)
   read(*,*) opt  

   call calGData(sF,sN,myGD)

   select case (opt)
   case (1) 
   sST = .true.
   do j = 1, sF
       print *
       print *, '-------------------------------'
       print *, ' Layer level:',j, ' sST=',sST
       print *, ' Layer Config.:', sN(1:sF)
       print *, ' Local Hosb Information:'       
       print *, " ----------------------------------------------------------------"
       print *, "  node  sN      pBlk     pDim     pSize     pStart   pEnd   pLen "
       print *, " -----------------------------------------------------------------"

       do i = 0, nNodes-1      
         call calGDMNode(i,nNodes,myGD,myNode)      
         call calMData(myGd, myNode, myData)
         call calHosbGrid(sST,myGD,myData,myHosb)
         write(*,100) i,sN(j),myData%pBlk(j),myData%pDim(j),myHosb%pSize(j), &
          myHosb%pStart(j),myHosb%pEnd(j),myHosb%pLen, myHosb%pMaxSize
       end do
   end do

   case (2)
   sST = .true.
   do j = 1, sF
       print *
       print *, '-------------------------------'
       print *, ' Layer level:',j, ' sST=',sST
       print *, ' Layer Config.:', sN(1:sF)
       print *, ' Global Hosb Information:'
       print *, " ----------------------------------------------------------------"
       print *, "  node  sN      gBlk     gDim     gSize     gStart   gEnd   gLen "
       print *, " -----------------------------------------------------------------"

       do i = 0, nNodes-1      
         call calGDMNode(i,nNodes,myGD,myNode)      
         call calMData(myGd, myNode, myData)
         call calHosbGrid(sST,myGD,myData,myHosb)
         write(*,100) i,sN(j),myGD%gBlk(j),myGD%gDim(j),myHosb%gSize(j), &
          myHosb%gStart(j),myHosb%gEnd(j),myHosb%gLen
       end do
   end do

   case (3)
   sST = .true.
   do j = 1, sF
       print *
       print *, '-------------------------------'
       print *, ' Layer level:',j, ' sST=',sST
       print *, ' Layer Config.:', sN(1:sF)
       print *, ' Global Hosb Position Information:'

       print *, " ----------------------------------------------------------------"
       print *, "  node  sN      gStart     gEnd   pStart     pEnd   gPos   gMaxSize "
       print *, " -----------------------------------------------------------------"

       do i = 0, nNodes-1      
         call calGDMNode(i,nNodes,myGD,myNode)      
         call calMData(myGd, myNode, myData)
         call calHosbGrid(sST,myGD,myData,myHosb)
         write(*,100) i,sN(j),myHosb%gStart(j),myHosb%gEnd(j),myHosb%pStart(j), &
          myHosb%pEnd(j),myHosb%gPos(j),myHosb%gMaxSize
       end do
   end do

   case (4) 
   sST = .false.
   do j = 1, sF
       print *
       print *, '-------------------------------'
       print *, ' Layer level:',j, ' sST=',sST
       print *, ' Layer Config.:', sN(1:sF)
       print *, ' Local Hosb Information:'       
       print *, " ----------------------------------------------------------------"
       print *, "  node  sN      pBlk     pDim     pSize     pStart   pEnd   pLen "
       print *, " -----------------------------------------------------------------"

       do i = 0, nNodes-1      
         call calGDMNode(i,nNodes,myGD,myNode)      
         call calMData(myGd, myNode, myData)
         call calHosbGrid(sST,myGD,myData,myHosb)
         write(*,100) i,sN(j),myData%pBlk(j),myData%pDim(j),myHosb%pSize(j), &
          myHosb%pStart(j),myHosb%pEnd(j),myHosb%pLen, myHosb%pMaxSize
       end do
   end do

   case (5)
   sST = .false.
   do j = 1, sF
       print *
       print *, '-------------------------------'
       print *, ' Layer level:',j, ' sST=',sST
       print *, ' Layer Config.:', sN(1:sF)
       print *, ' Global Hosb Information:'
       print *, " ----------------------------------------------------------------"
       print *, "  node  sN      gBlk     gDim     gSize     gStart   gEnd   gLen "
       print *, " -----------------------------------------------------------------"

       do i = 0, nNodes-1      
         call calGDMNode(i,nNodes,myGD,myNode)      
         call calMData(myGd, myNode, myData)
         call calHosbGrid(sst,myGD,myData,myHosb)
         write(*,100) i,sN(j),myGD%gBlk(j),myGD%gDim(j),myHosb%gSize(j), &
          myHosb%gStart(j),myHosb%gEnd(j),myHosb%gLen
       end do
   end do

   case (6)
   sST = .false.
   do j = 1, sF
       print *
       print *, '-------------------------------'
       print *, ' Layer level:',j, ' sST=',sST
       print *, ' Layer Config.:', sN(1:sF)
       print *, ' Global Hosb Position Information:'

       print *, " ----------------------------------------------------------------"
       print *, "  node  sN      gStart     gEnd   pStart     pEnd   gPos   gMaxSize "
       print *, " -----------------------------------------------------------------"

       do i = 0, nNodes-1      
         call calGDMNode(i,nNodes,myGD,myNode)      
         call calMData(myGd, myNode, myData)
         call calHosbGrid(sST,myGD,myData,myHosb)
         write(*,100) i,sN(j),myHosb%gStart(j),myHosb%gEnd(j),myHosb%pStart(j), &
          myHosb%pEnd(j),myHosb%gPos(j),myHosb%gMaxSize
       end do
   end do

   end select

   100 format (2(I5,1X),7(I8,1X))
   print *, '*********   Finish Testing  *******'

end program
!cccccccccccccccccccccccccccccccccccccccccc



