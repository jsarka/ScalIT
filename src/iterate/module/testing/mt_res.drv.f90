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

   integer :: nNodes, sF, sN(FMAX)
   Type(GDataInfo) :: myGd
   TYPE(MNodeInfo) :: myNode
   TYPE(MDataInfo) :: myData
   TYPE(MGridInfo) :: myRES

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
   do j = 1, sF
       print *
       print *, '-------------------------------'
       print *, ' Layer level:',j
       print *, ' Layer Config.:', sN(1:sF)
       print *, ' Local RES Information:'       
       print *, " ----------------------------------------------------------------"
       print *, "  node  sN      pBlk     pDim     pSize     pStart   pEnd   pLen "
       print *, " -----------------------------------------------------------------"

       do i = 0, nNodes-1      
         call calGDMNode(i,nNodes,myGD,myNode)      
         call calMData(myGd, myNode, myData)
         call calRESGrid(myGD,myData,myRES)
         write(*,100) i,sN(j),myData%pBlk(j),myData%pDim(j),myRES%pSize(j), &
          myRES%pStart(j),myRES%pEnd(j),myRES%pLen, myRES%pMaxSize
       end do
   end do

   case (2)
   do j = 1, sF
       print *
       print *, '-------------------------------'
       print *, ' Layer level:',j
       print *, ' Layer Config.:', sN(1:sF)
       print *, ' Global RES Information:'
       print *, " ----------------------------------------------------------------"
       print *, "  node  sN      gBlk     gDim     gSize     gStart   gEnd   gLen "
       print *, " -----------------------------------------------------------------"

       do i = 0, nNodes-1      
         call calGDMNode(i,nNodes,myGD,myNode)      
         call calMData(myGd, myNode, myData)
         call calRESGrid(myGD,myData,myRES)
         write(*,100) i,sN(j),myGD%gBlk(j),myGD%gDim(j),myRES%gSize(j), &
          myRES%gStart(j),myRES%gEnd(j),myRES%gLen
       end do
   end do

   case (3)
   do j = 1, sF
       print *
       print *, '-------------------------------'
       print *, ' Layer level:',j
       print *, ' Layer Config.:', sN(1:sF)
       print *, ' Global RES Position Information:'
       print *, " ----------------------------------------------------------------"
       print *, "  node  sN      gStart     gEnd   pStart     pEnd   gPos   gMaxSize "
       print *, " -----------------------------------------------------------------"

       do i = 0, nNodes-1      
         call calGDMNode(i,nNodes,myGD,myNode)      
         call calMData(myGd, myNode, myData)
         call calRESGrid(myGD,myData,myRES)
         write(*,100) i,sN(j),myRES%gStart(j),myRES%gEnd(j),myRES%pStart(j), &
          myRES%pEnd(j),myRES%gPos(j),myRES%gMaxSize
       end do
   end do
   end select

   100 format (2(I5,1X),7(I8,1X))
   print *, '*********   Finish Testing  *******'

end program
!cccccccccccccccccccccccccccccccccccccccccc



