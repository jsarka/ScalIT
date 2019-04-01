!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       IO Subroutines in MPI environment            c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadInitData()
    double precision :: t1, t2
    
    loadInitData = .false.

    if (id==rootID) then
         t1 = MPI_WTime()
         write(*,100), ' Loading H0 data from:', fH0
    end if
    if ( .NOT. loadH0())  return
    
    if ( sNDVR ) then
        if (id==rootID)   &
            write(*,100) ' Loading Out-most Layer data from:', fOUTH
        if ( .NOT. loadOutH()) return
        RES(1:myRES%pSize(sF))=0.0D0; ResSeq(1:myRES%pSize(1))=0.0D0;
        Eig0(1:pmax)=0.0D0
    else
        if (id==rootID)   &
           write(*,100) ' Loading RES data from:', fRES
        if (.NOT. loadRES()) return
    end if

    if (id==rootID)  then
        select case (sJOB)
        case (JOB_RES1,JOB_RES2)
          write(*,100) ' Loading absorption potential from:', fAPP

        case (JOB_CRP,JOB_CRP1,JOB_CRP2)
          write(*,100) ' Loading product absorption potential from:', fAPP
          write(*,100) ' Loading reactant absorption potential from:', fAPR

        case default
          print *, ' No absorption potential is needed for the calculation!'
        end select
    end if

    if (.NOT. loadAP()) return
    
    if (totalDep) then
        if (id==rootID)   &
           write(*,100) ' Loading Coord. Dep data from:', fDEP
        if (.NOT. loadDep()) return
    end if 

   if (id==rootID) then
      t2=MPI_WTime()
      write(*,*) "-------------------------------------------------"
      write(*,*) "  Time for loading all data (s): ", t2-t1
      write(*,*) "-------------------------------------------------"
      write(*,*) 
   end if

   loadInitData = .true.

   100 FORMAT(2x, A, A)
end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadH0()
   logical :: loadData, loadData_CX
   integer :: ierr
    
   if (sCX) then
      if (id==rootID)  loadH0 = loadData_CX(myH0%mLen,H0CX,srMode,fH0)
      call MPI_BCAST(H0CX,myH0%mLen,MPI_DOUBLE_COMPLEX,rootID,    &
                     MPI_COMM_WORLD,ierr)
   else
      if (id==rootID)  loadH0 = loadData(myH0%mLen,H0,srMode,fH0)
      call MPI_BCAST(H0,myH0%mLen,MPI_DOUBLE_PRECISION,rootID,    &
                     MPI_COMM_WORLD, ierr)
   end if

   loadH0 = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadRES()
   integer :: ierr

   call MLoadDataDiag(MPI_COMM_WORLD,fRES, myRES%gPos(sF),   &
               myconf%gDim(sF), myData%pDim(sF), sN(sF), RES, ierr)

   call MLoadData(MPI_COMM_WORLD,fRES,myRES%gPos(1),myRES%pSize(1), &
                  ResSeq,ierr)

   loadRES = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadEig()
   integer :: ierr

   call MLoadData(MPI_COMM_WORLD,fEig,myRES%gPos(1),myRES%pSize(1), &
                  Eig0,ierr)

   loadEig = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadRESEig()
   
   loadRESEig = loadRES()

   if (loadRESEig)  loadRESEig = loadEig()

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveEig()
   integer :: ierr

   call MSaveData(MPI_COMM_WORLD,fEig,myRES%gPos(1),myRES%pSize(1), &
                  EIG0, ierr)

   saveEig = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadAP()
     logical :: loadData

     select case (sJOB)
     case (JOB_RES1,JOB_RES2)
         loadAP = loadAP0()

     case (JOB_CRP,JOB_CRP1,JOB_CRP2) 
         loadAP = loadAPP()       
         if (loadAP) then
            loadAP = loadAPR()
            AP(1:myRES%pSize(1))=APP(1:myRES%pSize(1))+APR(1:myRES%pSize(1))
         end if

     case default
         loadAP = .true.
     end select

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadAP0()
   integer :: ierr

   call MLoadData(MPI_COMM_WORLD,fAPP,myRES%gPos(1),myRES%pSize(1), &
                  AP, ierr)

   loadAP0 = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadAPP()
   integer :: ierr

   call MLoadData(MPI_COMM_WORLD,fAPP,myRES%gPos(1),myRES%pSize(1), &
                  APP, ierr)

   loadAPP = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadAPR()
   integer :: ierr

   call MLoadData(MPI_COMM_WORLD,fAPR,myRES%gPos(1),myRES%pSize(1),  &
                  APR, ierr)

   loadAPR = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveHOSB()

   integer :: i, fh, ierr
   
   if (sST) then
      call MPI_FILE_OPEN(MPI_COMM_WORLD,fHOSB,MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                         MPI_INFO_NULL,fh,ierr)

      if (sCX) then
         do i=1,sF
            if (myNode%nodNum(i)>1) then
               call MWriteDataGrid_CX(fh,cxSize,myHOSB%gPos(i),myconf%gDim(i), &
                        myData%pDim(i),sN(i),HOSBCX(myHOSB%pStart(i)), ierr)
            else
               call MWriteData_CX(fh, cxSize,myHOSB%gPos(i),myHOSB%pSize(i),   &
                        HOSBCX(myHOSB%pStart(i)), ierr)
            end if
         end do
      else
         do i=1,sF
            if (myNode%nodNum(i)>1) then
               call MWriteDataGrid(fh, dbSize,myHOSB%gPos(i),myconf%gDim(i),   &
                         myData%pDim(i),sN(i),HOSB(myHOSB%pStart(i)), ierr)
            else
               call MWriteData(fh, dbSize,myHOSB%gPos(i),myHOSB%pSize(i),      &
                        HOSB(myHOSB%pStart(i)), ierr)
            end if
         end do
      end if

      call MPI_FILE_CLOSE(fh, ierr)

      saveHOSB = (ierr==0)

   else

      saveHOSB = .true.

   end if

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadHOSB()

   integer :: i, fh, ierr

   if (sST) then

      call MPI_FILE_OPEN(MPI_COMM_WORLD,fHOSB,MPI_MODE_RDONLY, &
                      MPI_INFO_NULL,fh,ierr)

      if (sCX) then
         do i=1,sF
            if (myNode%nodNum(i)>1) then
               call MReadDataGrid_CX(fh, cxSize,myHOSB%gPos(i),myconf%gDim(i),&
                      myData%pDim(i),sN(i),HOSBCX(myHOSB%pStart(i)), ierr)
            else
               call MReadData_CX(fh, cxSize,myHOSB%gPos(i),myHOSB%pSize(i),   &
                        HOSBCX(myHOSB%pStart(i)), ierr)
            end if
         end do
      else
         do i=1,sF
            if (myNode%nodNum(i)>1) then
               call MReadDataGrid(fh, dbSize,myHOSB%gPos(i),myconf%gDim(i),   &
                        myData%pDim(i),sN(i),HOSB(myHOSB%pStart(i)), ierr)
            else
               call MReadData(fh, dbSize,myHOSB%gPos(i),myHOSB%pSize(i),      &
                        HOSB(myHOSB%pStart(i)), ierr)
            end if
         end do
      end if

      call MPI_FILE_CLOSE(fh, ierr)

      loadHOSB  = (ierr==0)

   else

      loadHOSB = .true.
      
   end if

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveVOSB()

   integer :: i, fh, ierr

   call MPI_FILE_OPEN(MPI_COMM_WORLD,fVOSB,MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                      MPI_INFO_NULL,fh,ierr)

   do i=1,sF
      if (myNode%myID(i)==0)     &
         call MWriteData(fh, dbSize,myVOSB%gPos(i),myVOSB%pSize(i),  &
                        VOSB(myVOSB%pStart(i)), ierr)  
   end do  

   call MPI_FILE_CLOSE(fh, ierr)

   saveVOSB = (ierr==0)

end  function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadVOSB()  

   integer :: i, fh, ierr

   call MPI_FILE_OPEN(MPI_COMM_WORLD,fVOSB,MPI_MODE_RDONLY, &
                      MPI_INFO_NULL,fh,ierr)

   do i=1,sF
      if (myNode%myID(i)==0)     &
         call MReadData(fh, dbSize,myVOSB%gPos(i),myVOSB%pSize(i),  &
                        VOSB(myVOSB%pStart(i)), ierr)  
      if (myNode%nodNum(i)>1)   &
         call MPI_BCAST(VOSB(myVOSB%pStart(i)), myVOSB%pSize(i),    &
                        MPI_DOUBLE_PRECISION,0,myNode%commID(i), ierr)  
   end do  

   call MPI_FILE_CLOSE(fh, ierr)

   loadVOSB  = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveVOSBEig()

     saveVOSBEig = saveVOSB()

     if (saveVOSBEig) saveVOSBEig = saveEig()

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadVOSBEig()

     loadVOSBEig = loadVOSB()

     if (loadVOSBEig) loadVOSBEig = loadEig()

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadDep()
   integer :: i, fh, ierr

   loadDep = .true.
   if (.not. totalDep) return

   loadDep=.true.
   call MPI_FILE_OPEN(MPI_COMM_WORLD,fDEP,MPI_MODE_RDONLY,&
                      MPI_INFO_NULL,fh,ierr)

   if (sCX) then
      do i=1,sF
         if (sDep(i)) then
            if (myNode%nodNum(i)>1) then
               call MReadDataDiag_CX(fh, cxSize,myDep%gPos(i),myconf%gDim(i),&
                      sN(i), DEPCX(myDep%pStart(i)), ierr)
            else
               call MReadData_CX(fh, cxSize,myDep%gPos(i),myDep%pSize(i),&
                        DEPCX(myDEP%pStart(i)), ierr)
            end if
         end if
      end do
   else
      do i=1,sF
         if (sDep(i)) then
            if (myNode%nodNum(i)>1) then
               call MReadDataDiag(fh, dbSize,myDep%gPos(i),myconf%gDim(i),  &
                        sN(i), DEP(myDEP%pStart(i)), ierr)
            else
               call MReadData(fh, dbSize,myDEP%gPos(i),myDEP%pSize(i),  &
                        DEP(myDEP%pStart(i)), ierr)
            end if
         end if
      end do
   end if

   call MPI_FILE_CLOSE(fh, ierr)

   loadDep  = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadOutH()
   integer :: ierr

   if (sCX) then
       call MLoadDataGrid_CX(MPI_COMM_WORLD,fOUTH,myRES%gPos(sF),  &
              myconf%gDim(sF),myData%pDim(sF),sN(sF),OUTHCX, ierr)
   else
       call MLoadDataGrid(MPI_COMM_WORLD,fOutH,myRES%gPos(sF),     &
              myconf%gDim(sF), myData%pDim(sF),sN(sF),OUTH, ierr)
   end if

   loadOutH = (ierr==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


