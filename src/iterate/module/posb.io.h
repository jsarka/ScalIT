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
          if (.NOT. loadAP()) return

        case (JOB_CRP,JOB_CRP1,JOB_CRP2)
          write(*,100) ' Loading product absorption potential from:', fAPP
          write(*,100) ' Loading reactant absorption potential from:', fAPR
          if (loadAPP()) then
             if (loadAPR()) then
                AP(1:plen(1))=APP(1:plen(1))+APR(1:plen(1))
             else
                return
             end if
          else
             return
          end if

        case default
          print *, ' No absorption potential is needed for the calculation!'
        end select
    end if

    if (totalDep) then
        if (id==rootID) &
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
logical function loadRESEig()
   
   loadRESEig = loadRES()

   if (loadRESEig)  loadRESEig = loadEig()

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

