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

   integer :: i, j

   logical :: sNDVR
   integer :: nNodes, sF, sN(FMAX)
   Type(GDataInfo) :: myGd
   TYPE(MNodeInfo) :: myNode
   TYPE(MDataInfo) :: myData

   print *,'  =================================='
   print *,'        Testing MOSBType:calMData   '
   print *, ' =================================='
   print *
   print *, ' Read input parameters from standard input. '

   read(*,*) nNodes, sF
   read(*,*) sN(1:sF)
  
   call calGData(sF,sN,myGD)

   do j = 1, sF
       print *
       print *, '-------------------------------'
       print *, ' Layer level:',j
       print *, ' Layer Config.:', sN(1:sF)
       print *, ' Block Information'
       print *, " ------------------------------------------------------"
       print *, "     node        Start       End    pSize    pBlk*pDim "
       print *, " ------------------------------------------------------"

       do i = 0, nNodes-1      
         call calGDMNode(i,nNodes,myGD,myNode)      
         call calMData(myGd, myNode, myData)
         write(*,100) i,myData%gBStart(j),myData%gBEnd(j),myData%pBlk(j),myData%pLen(j)
      end do
   end do

   do j = 1, sF
       print *
       print *, '-------------------------------'
       print *, ' Layer level:',j
       print *, ' Layer Config.:', sN(1:sF)
       print *, ' Dimension Information'
       print *, " ------------------------------------------------------"
       print *, "     node        Start       End    pSize    pMaxLen "
       print *, " ------------------------------------------------------"

       do i = 0, nNodes-1      
         call calGDMNode(i,nNodes,myGD,myNode)      
         call calMData(myGd, myNode, myData)
         write(*,100) i,myData%gDStart(j),myData%gDEnd(j),myData%pDim(j),myData%pMaxLen
      end do
   end do

   print *, '*********   Finish Testing  *******'

   100 FORMAT(2X,I6, 1X, 2I12, I6, 2I12)

end program
!cccccccccccccccccccccccccccccccccccccccccc



