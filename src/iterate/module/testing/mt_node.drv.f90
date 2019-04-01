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

   print *,'  =================================='
   print *,'        Testing MOSBType:calMNode   '
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
       print *, " -------------------------------------------------------------------"
       print *, "     node     myID   GrpID     commID       #Grp      Num    nStart "
       print *, " -------------------------------------------------------------------"
       do i = 0, nNodes-1
         call calGDMNode(i,nNodes,myGd,myNode)
         write(*,100) i, myNode%myID(j),myNode%grpID(j),myNode%commID(j),   &
                        myNode%nGroup(j),myNode%nodNum(j),myNode%nodIDStart(j)
      end do
   end do
   print *
   print *, 'spevel:', myNode%splevel
   print *, 'Balance Flag:', myNode%lbFlag
   print *, '*********   Finish Testing  *******'
   100 FORMAT(2X,3(I7,1X),I14, 1X, 3(I7,1X))

end program
!cccccccccccccccccccccccccccccccccccccccccc



