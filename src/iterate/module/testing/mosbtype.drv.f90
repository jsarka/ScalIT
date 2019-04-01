!ccccccccccccccccccccccccccccccccccccccccccccccc
!c        Testing program for mosbtype         c
!ccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccc
!c   comment out MPI_COMM_SPLIT function call   c
!c    at line 49 in file mosbtype.cal.h         c
!cccccccccccccccccccccccccccccccccccccccccccccccc
program mosbType_test    
   use mosbtype
   implicit none

   integer :: opt, i

   logical :: sST, sNDVR
   TYPE(CConvSimple) :: mySConv
   TYPE(CConv)     :: myConv
   TYPE(COsbw)     :: myOsbw
   Type(BConfInfo) :: myconf
   TYPE(MNodeInfo) :: myNode
   Type(MDataInfo) :: myData
   TYPE(CSeqInfo)  :: myH0
   TYPE(MGridInfo) :: myHOSB, myVOSB, myRES, myDep
   logical :: sDep(FMAX)

   print *,'  =================================='
   print *,'        Testing MOSBType            '
   print *, ' =================================='
   print *
   print *, ' Read input parameters from standard input. '

   call readInput(opt, myconf, mysconv, myconv, myosbw, sDep)
  
   select case(opt)
   case (1)
      print *
      print *, ' Print Simple Conv'
      call printSConv(mysconv)

      print *
      print *, ' Print Conv'
      call printCConv(myconv)

      print *
      print *, ' Print OSBW'
      call printOSBW(myosbw)
 
      print *
      print *, ' Print Basic Config'
      call printBInfo(myconf)

   case (2)
      print *
      call printBInfo(myconf)
      do i = 1, myconf%nNodes
         myconf%id=i-1
         call calMNode(myconf, mynode)
         print *, ' Print MyNode Info'
         call printMNode(myconf, myNode)
       end do
   

   case (3)
      print *
      call printBInfo(myconf)
      do i = 1, myconf%nNodes
         myconf%id=i-1
         call calMNode(myconf, mynode)
         call calMData(myconf,myNode, myData)
         print *, ' Print MyData Info'
         call printMData(myconf, myData)
      end do

   case (4)
      print *
      call printBInfo(myconf)
      do i = 1, myconf%nNodes
         myconf%id=i-1
         call calH0Seq(.true.,myconf, myH0)         
         print *, ' Print MyH0 Info: TRUE'
         call printSeqInfo(myH0)
         call calH0Seq(.false.,myconf, myH0)
         print *, ' Print MyH0 Info: FALSE'
         call printSeqInfo(myH0)
       end do

   case (5)
      print *
      call printBInfo(myconf)
      do i = 1, myconf%nNodes
         myconf%id=i-1
         call calMNode(myconf, mynode)
         call calMData(myconf,myNode, myData)
         call calHOSBGrid(ssT,myconf,mynode,myData,myHOSB)
         print *, ' Print MyHOSB Info'
         call printMGrid(myconf,myHOSB)
       end do

   case (6)
      print *
      call printBInfo(myconf)
      do i = 1, myconf%nNodes
         myconf%id=i-1
         call calMNode(myconf, mynode)
         call calMData(myconf,myNode, myData)
         call calVOSBGrid(myconf,mynode,myData,myVOSB)
         print *, ' Print MyVOSB Info'
         call printMGrid(myconf,myVOSB)
       end do

   case (7)
      print *
      call printBInfo(myconf)
      do i = 1, myconf%nNodes
         myconf%id=i-1
         call calMNode(myconf,mynode)
         call calMData(myconf,myNode, myData)
         call calRESGrid(myconf,mynode,myData,myRES)
         print *, ' Print MyRES Info'
         call printMGrid(myconf,myRES)
      end do

   case (8)
      print *
      call printBInfo(myconf)
      do i = 1, myconf%nNodes
         myconf%id=i-1
         call calMNode(myconf,mynode)
         call calMData(myconf,myNode,myData)
         call calDepGrid(myconf,mynode,myData,sDep,myDep)
         print *, ' Print MyDep Info'
         call printMGrid(myconf,myDep)
      end do

   end select

end program
!cccccccccccccccccccccccccccccccccccccccccc

subroutine readInput(opt, myconf, mysconv, myconv, myosbw, sDep )
   use mosbtype
   implicit none
!   use mosbtype
   integer,           intent(OUT) :: opt
   TYPE(CConvSimple), intent(OUT) :: mySConv
   TYPE(CConv),       intent(OUT) :: myConv
   TYPE(COsbw),       intent(OUT) :: myOsbw
   Type(BConfInfo),   intent(OUT) :: myconf
   logical,           intent(OUT) :: sDep(FMAX)

   integer :: i

   read(*,*) opt
   read(*,*) mySconv
   read(*,*) myConv
   read(*,*) myosbw
   read(*,*) myconf%id, myconf%nNodes,myconf%sF
   read(*,*) (myconf%sN(i),i=1,myconf%sF)
   read(*,*) (sDep(i),i=1,myconf%sF)

end 


