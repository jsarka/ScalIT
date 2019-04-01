!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       IO Subroutines in MPI environment            c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c  Works without the parallel IO,
!c

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadInitDataSeq()
    
    loadInitDataSeq = .false.

    if (id==rootID)   &
         write(*,100), ' Loading H0 data from:', fH0
    if ( .NOT. loadH0())  return
    
    if ( sNDVR ) then
        if (id==rootID)   &
            write(*,100) ' Loading Out-most Layer data from:', fOUTH
        if ( .NOT. loadOutHSeq(OutH)) then
           print *, ' Error in loadOutHSeq !'
           return      
        end if
        RES(1:plen(sF))=0.0D0;  ResSeq(1:plen(1))=0.0D0
        Eig0(1:pmax)=0.0D0
    else
        if (id==rootID)   &
           write(*,100) ' Loading RES and Eig0 data from:', fRES
        if (.NOT. loadRESSeq()) return
    end if

    if (id==rootID)  then
       select case (sJOB)
       case (JOB_RES1,JOB_RES2)
          write(*,100) ' Loading absorption potential from:', fAPP

       case (JOB_CRP,JOB_CRP1,JOB_CRP2)
          write(*,100) ' Loading product absorption potential from:', fAPP
          write(*,100) ' Loading reactant absorption potential from:', fAPR
       end select
    end if
    if (.NOT. loadAPSeq()) return
   
    if (totalDep) then
        if (id==rootID)   &
           write(*,100) ' Loading Coord. Dep data from:', fDEP
        if (.NOT. loadDepSeq()) return
    end if 

   loadInitDataSeq = .true.

   100 FORMAT(2x, A, A)
end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadRESSeq()

   loadRESSeq = loadDiag(srMode,fRes, RES)
   if (loadRESSeq) loadREsSeq = loadSeq(srMode,fRes,RESSeq)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadAPSeq()
     logical :: loadData

     select case (sJOB)
     case (JOB_RES1, JOB_RES2)
         loadAPSeq = loadAPSeq0()

     case (JOB_CRP,JOB_CRP1,JOB_CRP2)
         loadAPSeq = loadAPPSeq()
         if (loadAPSeq) then
            loadAPSeq = loadAPRSeq()
            AP(1:plen(1))=APP(1:plen(1))+APR(1:plen(1))
         else
            return
         end if

     case default
         loadAPSeq = .true.
     end select

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadAPSeq0()

   loadAPSeq0 = loadSeq(srMode,fAPP, AP)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadAPRSeq()

   loadAPRSeq = loadSeq(srMode,fAPR, APR)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadAPPSeq()

   loadAPPSeq = loadSeq(srMode,fAPP, APP)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


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
logical function loadOutHSeq(locH0)
   double precision, intent(OUT) :: locH0(locDim(id+1,SF),SN(SF),SN(SF))

   integer, parameter :: OUTH_TAG = 3000
   double precision, allocatable :: gOutH(:)
   integer :: len(nNodes),disp(nNodes)
   integer :: glen, loclen, recLen, ierr, info, recInd    
   integer :: i, j, k

   loadOutHSeq=.false.  
   glen=myconf%gDim(sF); locLen=locDim(id+1,sF)
   if (rootID==id) then
      allocate(gOutH(glen))
   else
      allocate(gOutH(1))
   end if

   inquire(IOLENGTH=recLen) gOutH
   if (id==rootID)  &
     open(99,file=fOutH,STATUS='OLD',FORM='UNFORMATTED',  &
          ACCESS='direct',RECL=recLen,IOSTAT=info)

   call MPI_BCast(info,1,MPI_INTEGER,rootID, MPI_COMM_WORLD, ierr)

   if (info==0) then 
      do i = 1, nNodes
	 len(i)=locDim(i,sF)
      end do
      disp(1)=0
      do i=1,nNodes-1
	  disp(i+1)=disp(i)+len(i)
      end do
      if ((disp(nNodes)+len(nNodes))==glen) then
        recInd=1       
	do i = 1, SN(SF)
           do j = 1, sN(SF)
	      if (id==rootID)  &
                  read(99,rec=recInd) gOutH
	      call MPI_SCATTERV(gOutH,len,disp,MPI_DOUBLE_PRECISION,locH0(1,j,i), &
                    locLen,MPI_DOUBLE_PRECISION, rootID,MPI_COMM_WORLD,ierr  ) 
              recInd=recInd+1
           end do
        end do
        loadOutHSeq=.true.
      else
	if (id==rootID) &  
	  print *, ' Error for data size in load OutH data'
      end if
      if (id==rootID)  close(99)      
   end if

   if (allocated(gOutH))  deallocate(gOutH)

   loadOutHSeq = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadSeq(rdMode, fname, dat1)
   logical, intent(IN) :: rdMode
   character(len=*), intent(IN) :: fname
   double precision, intent(OUT):: dat1(myRES%pSize(1))

   integer, parameter :: GDAT_TAG = 1000
   integer :: ierr, totlen, ind,status(MPI_STATUS_SIZE)
   double precision, allocatable :: gDat(:)

   integer :: i

   totlen=myRes%gLen

   if (id==rootID) then
      allocate(gDat(totlen))
      call loadDataDir(totLen,gDat,fname)
   end if

   if (id==rootID) then
      do i = 0, nNodes-1
         if (i==0) then
            ind = 1
         else
	     ind = ind + bNum(i,1)*locDim(i,1)*sN(1)
         end if
         if (i/=rootid) then
            call MPI_Send(gDat(ind),bNum(i+1,1)*locDim(i+1,1)*sN(1),  &
                 MPI_DOUBLE_PRECISION, i, GDAT_TAG, MPI_COMM_WORLD, ierr)
         else
	    dat1(1:myRES%pSize(1))=gDat(ind:ind+myRES%pSize(1)-1)
         end if
      end do
   else
      call MPI_Recv(dat1,myRES%pSize(1),MPI_DOUBLE_PRECISION,  &
                 rootid,GDAT_TAG,MPI_COMM_WORLD,status,ierr)
   end if

   if (allocated(gDat))  deallocate(gDat)

   loadSeq = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveSeq(svMode,fname,dat1)
   logical, intent(IN) :: svMode
   character(len=*), intent(IN) :: fname
   double precision, intent(IN):: dat1(myRES%pSize(1))

   integer, parameter :: GDAT_TAG = 1000
   integer :: ierr, totlen, ind,status(MPI_STATUS_SIZE)
   double precision, allocatable :: gDat(:)

   integer :: i

   totlen=myRES%gLen

   if (id==rootID)   allocate(gDat(totlen))

   if (id==rootID) then
      do i = 0, nNodes-1
         if (i==0) then
            ind = 1
         else
 	    ind = ind + bNum(i,1)*locDim(i,1)*sN(1)
         end if
         
         if (i/=rootID) then
            call MPI_Recv(gDat(ind),bNum(i+1,1)*locDim(i+1,1)*sN(1),  &
                 MPI_DOUBLE_PRECISION,i,GDAT_TAG,MPI_COMM_WORLD,status,ierr)
         else
	    gDat(ind:ind+myRES%pSize(1)-1)=dat1(1:myRES%pSize(1))
         end if
      end do
   else
      call MPI_Send(dat1,myRES%pSize(1),MPI_DOUBLE_PRECISION,  &
                 rootid, GDAT_TAG, MPI_COMM_WORLD, ierr)
   end if

   if (id==rootID) then
      call saveDataDir(totLen,gDat,svMode,fname)
      deallocate(gDat)
   end if

   saveSeq = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadDiag(rdMode, fname, dat1)
   logical, intent(IN) :: rdMode
   character(len=*), intent(IN)  :: fname
   double precision, intent(OUT) :: dat1(locDim(id+1,sF),sN(sF))

   integer, parameter :: DIAG_TAG = 2000
   integer :: ierr, totlen, ind,status(MPI_STATUS_SIZE)
   double precision, allocatable :: gDat(:,:)
   
   logical :: ltemp, loadData

   integer :: j, k
   totlen=myconf%gDim(sF)

   if (id==rootID) then
      allocate(gDat(totlen,sN(sF)))
!      call loadDataDir(totLen*sN(sF),gDat,rdMode,fname)
      ltemp = loadData(totLen*sN(sF),gDat,rdMode,fname)
      print *, 'load data from ', fname
      print *, 'data:',gdat
      print *
   end if

   do j = 1, sN(sF)
      if (id==rootID) then
         do k = 0, nNodes-1
            if (k==0) then
               ind = 1
            else
               ind = ind + locDim(k,sF)
            end if

            if (k/=rootID) then
               call MPI_Send(gDat(ind,j),locDim(k+1,sF),  &
                   MPI_DOUBLE_PRECISION,k,DIAG_TAG,MPI_COMM_WORLD,ierr)
            else
   	       dat1(1:locDim(id+1,sF),j)=gDat(ind:ind+locDim(id+1,sF)-1,j)
            end if
         end do
      else
         call MPI_Recv(dat1(1, j),locDim(id+1,sF),MPI_DOUBLE_PRECISION, &
		     rootid, DIAG_TAG, MPI_COMM_WORLD,status,ierr)
      end if

   end do

   if (id==rootID)  deallocate(gDat)

   loadDiag = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadDiag_CX(rdMode, fname, dat1)
   logical, intent(IN) :: rdMode
   character(len=*), intent(IN)  :: fname
   double complex, intent(OUT) :: dat1(locDim(id+1,sF),sN(sF))

   integer, parameter :: DIAG_TAG = 3000
   integer :: ierr, totlen, ind,status(MPI_STATUS_SIZE)
   double complex, allocatable :: gDat(:,:)
   
   logical :: ltemp, loadData_CX

   integer :: j, k
   totlen=myconf%gDim(sF)

   if (id==rootID) then
      allocate(gDat(totlen,sN(sF)))
!      call loadDataDir_CX(totLen*sN(sF),gDat,fname)
      ltemp = loadData_CX(totLen*sN(sF),gDat,rdMode,fname)
   end if

   do j = 1, sN(sF)
      if (id==rootID) then
         do k = 0, nNodes-1
            if (k==0) then
               ind = 1
            else
               ind = ind + locDim(k,sF)
            end if
            
            if (k/=rootID) then
               call MPI_Send(gDat(ind,j),locDim(k+1,sF),  &
                   MPI_DOUBLE_COMPLEX,k,DIAG_TAG,MPI_COMM_WORLD,ierr)
            else
   	       dat1(1:locDim(k+1,sF),j)=gDat(ind:ind+locDim(k+1,sF)-1,j)
            end if
         end do
      else
         call MPI_Recv(dat1(1, j),locDim(id+1,sF),MPI_DOUBLE_COMPLEX, &
		     rootid, DIAG_TAG, MPI_COMM_WORLD,status,ierr)
      end if

   end do

   if (id==rootID)  deallocate(gDat)

   loadDiag_CX = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveDiag(svMode, fname, dat1)
   logical, intent(IN) :: svMode
   character(len=*), intent(IN)  :: fname
   double precision, intent(IN) :: dat1(locDim(id+1,sF),sN(sF))

   integer, parameter :: DIAG_TAG = 4000
   integer :: ierr, totlen, ind,status(MPI_STATUS_SIZE)
   double precision, allocatable :: gDat(:,:)

   integer :: j, k
   totlen=myconf%gDim(sF)

   if (id==rootID) then
      allocate(gDat(totlen,sN(sF)))
   end if

   do j = 1, sN(sF)
      if (id==rootID) then
         do k = 0, nNodes-1
            if (k==0) then
               ind = 1
            else
               ind = ind + locDim(k,sF)
            end if
            
            if (k/=rootID) then
               call MPI_Recv(gDat(ind,j),locDim(k+1,sF),MPI_DOUBLE_PRECISION,&
                   k,DIAG_TAG,MPI_COMM_WORLD,status,ierr)
            else
   	       gDat(ind:ind+locDim(k+1,sF)-1,j) = dat1(1:locDim(k+1,sF),j)
            end if
         end do
      else
         call MPI_SEND(dat1(1, j),locDim(id+1,sF),MPI_DOUBLE_PRECISION, &
		     rootid, DIAG_TAG, MPI_COMM_WORLD,ierr)
      end if

   end do

   if (id==rootid) then
!       call saveDataDir(totLen*sN(sF),gDat,svMode,fname)
       call saveData(totLen*sN(sF),gDat,svMode,fname)
       print *, ' Save data in file:',fname
       print *, 'Save data:', gDat
       print *
       deallocate(gDat)
   end if

   saveDiag = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveDiag_CX(svMode, fname, dat1)
   logical, intent(IN) :: svMode
   character(len=*), intent(IN)  :: fname
   double complex, intent(IN) :: dat1(locDim(id+1,sF),sN(sF))

   integer, parameter :: DIAG_TAG = 5000
   integer :: ierr, totlen, ind,status(MPI_STATUS_SIZE)
   double complex, allocatable :: gDat(:,:)

   integer :: j, k
   totlen=myconf%gDim(sF)

   if (id==rootID) then
      allocate(gDat(totlen,sN(sF)))
   end if

   do j = 1, sN(sF)
      if (id==rootID) then
         do k = 0, nNodes-1
            if (k==0) then
               ind = 1
            else
               ind = ind + locDim(k,sF)
            end if
            
            if (k/=rootID) then
               call MPI_Recv(gDat(ind,j),locDim(k+1,sF),MPI_DOUBLE_COMPLEX,&
                   k,DIAG_TAG,MPI_COMM_WORLD,status,ierr)
            else
   	       gDat(ind:ind+locDim(k+1,sF)-1,j) = dat1(1:locDim(k+1,sF),j)
            end if
         end do
      else
         call MPI_SEND(dat1(1, j),locDim(id+1,sF),MPI_DOUBLE_COMPLEX, &
		     rootid, DIAG_TAG, MPI_COMM_WORLD,ierr)
      end if

   end do

   if (id==rootid) then
!       call saveDataDir_CX(totLen*sN(sF),gDat,svMode,fname)
       call saveData_CX(totLen*sN(sF),gDat,svMode,fname)
       deallocate(gDat)
   end if

   saveDiag_CX = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
