!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       IO Subroutines in MPI environment without parallel IO    c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadDepSeq()
   integer :: i, fh, ierr

   if (.not. totalDep) return

   call MPI_FILE_OPEN(MPI_COMM_WORLD,fDep,MPI_MODE_RDONLY,     &
                      MPI_INFO_NULL,fh,ierr)

   if (sCX) then
      do i=1,sF
         if (myNode%nodNum(i)>1) then
            call MReadDataDiag_CX(fh, cxSize,myDep%gPos(i),myconf%gDim(i),  &
                      sN(i), DEPCX(myDep%pStart(i)), ierr)
         else
            call MReadData_CX(fh, cxSize,myDep%gPos(i),myDep%pSize(i),&
                        DEPCX(myDEP%pStart(i)), ierr)
         end if
      end do
   else
      do i=1,sF
         if (myNode%nodNum(i)>1) then
            call MReadDataDiag(fh, dbSize,myDep%gPos(i),myconf%gDim(i),sN(i),&
                        DEP(myDEP%pStart(i)), ierr)
         else
            call MReadData(fh, dbSize,myDEP%gPos(i),myDEP%pSize(i),  &
                        DEP(myDEP%pStart(i)), ierr)
         end if
      end do
   end if

   call MPI_FILE_CLOSE(fh, ierr)

   loadDepSeq = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadSeq(layer, fname, dat1)
   integer, intent(IN) :: layer
   character(len=*), intent(IN) :: fname
   double precision, intent(OUT):: dat1(myRES%pSize(layer))

   integer :: ierr

   call PLoadDataSeq(id,rootID,nNodes,locDim(1:nNodes,layer), MPI_COMM_WORLD,   &
           fname,myRES%pSize(layer),dat1,ierr)
   
   loadSeq = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadSeq_CX(layer,fname, dat1)
   integer, intent(IN) :: layer
   character(len=*), intent(IN) :: fname
   double complex, intent(OUT):: dat1(myRES%pSize(layer))

   integer :: ierr

   call PLoadDataSeq_CX(id,rootID,nNodes, locDim(1:nNodes,layer), MPI_COMM_WORLD,   &
           fname,myRES%pSize(layer),dat1,ierr)

   loadSeq_CX = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveSeq(layer,fname,dat1)
   integer, intent(IN) :: layer
   character(len=*), intent(IN) :: fname
   double precision, intent(IN):: dat1(myRES%pSize(layer))

   integer :: ierr

   call PSaveDataSeq(id,rootID,nNodes, locDim(1:nNodes,layer), MPI_COMM_WORLD,  &
           fname,myRES%pSize(layer),dat1,ierr)

   saveSeq = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveSeq_CX(layer,fname,dat1)
   integer, intent(IN) :: layer
   character(len=*), intent(IN) :: fname
   double complex, intent(IN):: dat1(myRES%pSize(layer))

   integer :: ierr

   call PSaveDataSeq_CX(id,rootID,nNodes, locDim(1:nNodes,layer), MPI_COMM_WORLD,  &
           fname,myRES%pSize(layer),dat1,ierr)

   saveSeq_CX = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadDiag(layer,fname, dat1)
   integer, intent(IN) :: layer
   character(len=*), intent(IN)  :: fname
   double precision, intent(OUT) :: dat1(locDim(id+1,layer),sN(layer))

   integer :: ierr
   
   call PLoadDataDiag(id,rootID,nNodes,locDim(1:nNodes,layer),  &
            MPI_COMM_WORLD,fname,locDim(id+1,layer),sN(layer),dat1,ierr)

   loadDiag = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadDiag_CX(layer, fname, dat1)
   integer, intent(IN) :: layer
   character(len=*), intent(IN)  :: fname
   double complex, intent(OUT) :: dat1(locDim(id+1,layer),sN(layer))

   integer :: ierr
   
   call PLoadDataDiag_CX(id,rootID,nNodes,locDim(1:nNodes,layer),  &
            MPI_COMM_WORLD,fname,locDim(id+1,layer),sN(layer),dat1,ierr)

   loadDiag_CX = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveDiag(layer,fname, dat1)
   integer, intent(IN) :: layer
   character(len=*), intent(IN)  :: fname
   double precision, intent(IN) :: dat1(locDim(id+1,layer),sN(layer))

   integer :: ierr

   call PSaveDataDiag(id,rootID,nNodes,locDim(1:nNodes,layer),  &
            MPI_COMM_WORLD,fname,locDim(id+1,layer),sN(layer),dat1,ierr)
 
   saveDiag = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveDiag_CX(layer,fname, dat1)
   integer, intent(IN) :: layer
   character(len=*), intent(IN)  :: fname
   double complex, intent(IN) :: dat1(locDim(id+1,layer),sN(layer))

   integer :: ierr

   call PSaveDataDiag_CX(id,rootID,nNodes,locDim(1:nNodes,layer),  &
            MPI_COMM_WORLD,fname,locDim(id+1,layer),sN(layer),dat1,ierr)
 
   saveDiag_CX = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadGrid(layer,fname, dat1)
   integer, intent(IN) :: layer
   character(len=*), intent(IN)  :: fname
   double precision, intent(OUT) :: dat1(locDim(id+1,layer),sN(layer),SN(layer))

   integer :: ierr

   call PLoadDataGrid(id,rootID,nNodes,locDim(1:nNodes,layer),  &
            MPI_COMM_WORLD,fname,locDim(id+1,layer),sN(layer),dat1,ierr)

   loadGrid = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadGrid_CX(layer, fname, dat1)
   integer, intent(IN) :: layer
   character(len=*), intent(IN)  :: fname
   double complex, intent(OUT) :: dat1(locDim(id+1,layer),sN(layer),SN(layer))

   integer :: ierr

   call PLoadDataGrid_CX(id,rootID,nNodes,locDim(1:nNodes,layer),  &
            MPI_COMM_WORLD,fname,locDim(id+1,layer),sN(layer),dat1,ierr)

   loadGrid_CX = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveGrid(layer,fname, dat1)
   integer, intent(IN) :: layer
   character(len=*), intent(IN)  :: fname
   double precision, intent(IN) :: dat1(locDim(id+1,layer),sN(layer),sN(layer))

   integer :: ierr

   call PSaveDataGrid(id,rootID,nNodes,locDim(1:nNodes,layer),  &
            MPI_COMM_WORLD,fname,locDim(id+1,layer),sN(layer),dat1,ierr)

   saveGrid = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveGrid_CX(layer,fname, dat1)
   integer, intent(IN) :: layer
   character(len=*), intent(IN)  :: fname
   double complex, intent(IN) :: dat1(locDim(id+1,layer),sN(layer),sN(layer))

   integer :: ierr

   call PSaveDataGrid_CX(id,rootID,nNodes,locDim(1:nNodes,layer),  &
            MPI_COMM_WORLD,fname,locDim(id+1,layer),sN(layer),dat1,ierr)

   saveGrid_CX = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

