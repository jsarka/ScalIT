!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadInitData()
    
    loadInitData = .false.

    write(*,100), ' Loading H0 data from:', H0File
    if ( .NOT. loadH0())  return
    
    if ( sNDVR ) then
        write(*,100) ' Loading Out-most Layer data from:', OUTHFile
        if ( .NOT. loadOutH()) return
        RES(1:myLen)=0.0D0; Eig0(1:myLen)=0.0D0
    else
        write(*,100) ' Loading RES and Eig0 data from:', RESFile
        if (.NOT. loadRESEig()) return
    end if

    select case (sJOB)

    case (JOB_RES1,JOB_RES2)   ! load APP, APR=0, AP=APP
       write(*,100) ' Loading Absorption Potential from:', APPFile
       if (.NOT. loadAPP()) return
       APR(1:mylen)=0.0D0; AP(1:mylen)=APP(1:mylen)

    case (JOB_CRP)         ! load APP, APR
       write(*,100) ' Loading product absorption potential from:', APPFile
       write(*,100) ' Loading reactant absorption potential from:', APRFile
       if (.NOT. loadAP()) return

    case default
       print *, '  No Absorption potential is used!'
    end select

    if (totalDep) then
        write(*,100) ' Loading Coord. Dep data from:', DEPFile
        if (.NOT. loadDep()) return
    end if 

   loadInitData = .true.

   100 FORMAT(2x, A, A)
end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadH0()
     logical :: loadData, loadData_CX

     if (sCX) then
        loadH0 = loadData_CX(myH0%mLen,H0CX,srMode,H0FILE)
     else
        loadH0 = loadData(myH0%mLen,H0,srMode,H0FILE)
     end if

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadRES()
     logical :: loadDataDir
     
     loadRES = loadDataDir(myLen,RES,RESFILE)
end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadEig()
     logical :: loadDataDir
     
     loadEIG = loadDataDir(myLen,Eig0,EIGFILE)
end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadRESEig()
     logical :: loadDataDir
     
     loadRESEig = loadDataDir(myLen,RES,RESFILE)
     EIG0(1:myLen) = RES(1:myLen)
end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveEig()
     logical :: saveDataDir
     
     saveEig = saveDataDir(myLen,Eig0,EIGFILE)
end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadAP()
     logical :: loadDataDir

     loadAP = loadDataDir(myLen,APP,APPFILE)

     if (loadAP) then
        loadAP = loadDataDir(myLen,APR,APRFILE)
        if (loadAP)  AP(1:myLen)=APP(1:myLen)+APR(1:myLen)
     end if

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadAPP()
     logical :: loadDataDir
     
     loadAPP = loadDataDir(myLen,APP,APPFILE)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadAPR()
     logical :: loadDataDir
    
     loadAPR = loadDataDir(myLen,APR,APRFILE)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadHOSB()
     logical :: loadData, loadData_CX, loadDataDir, loadDataDir_CX

     loadHOSB=.false.
     
     if ((sST) .AND. (sHOSB<0)) then     
        if (sHOSB<-SEQDIRNUM) then
           if (sCX) then
               loadHOSB = loadData_CX(myHOSB%mLen,HOSBCX,srMode,HOSBFILE)
           else
               loadHOSB = loadData(myHOSB%mLen,HOSB,srMode,HOSBFILE)
           end if
        else
           if (sCX) then
               loadHOSB = loadDataDir_CX(myHOSB%mLen,HOSBCX,HOSBFILE)
           else
               loadHOSB = loadDataDir(myHOSB%mLen,HOSB,HOSBFILE)
           end if
        end if
     end if

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveHOSB()
     logical :: saveData, saveData_CX, saveDataDir, saveDataDir_CX

     saveHOSB = .false.
     if ((sST).AND.(sHOSB>0)) then
        if (sHOSB>SEQDIRNUM) then
           if (sCX) then
              saveHOSB = saveData_CX(myHOSB%mLen,HOSBCX,srMode,HOSBFILE)
           else
              saveHOSB = saveData(myHOSB%mLen,HOSB,srMode,HOSBFILE)
           end if
        else
           if (sCX) then
              saveHOSB = saveDataDir_CX(myHOSB%mLen,HOSBCX,HOSBFILE)
           else
              saveHOSB = saveDataDir(myHOSB%mLen,HOSB,HOSBFILE)
           end if
        end if
     end if
end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadVOSBEig()
     logical :: loadData,loadDataDir

     loadVOSBEig = .false.

     if (sVOSB>=0) return

     if (sVOSB<-SEQDIRNUM) then
        loadVOSBEig = loadData(myVOSB%mLen,VOSB,srMode,VOSBFILE)
        if (loadVOSBEig) loadVOSBEig = loadData(myLen,Eig0,srMode,EIGFILE)

     else
        loadVOSBEig = loadDataDir(myVOSB%mLen,VOSB,VOSBFILE)
        if(loadVOSBEig) loadVOSBEig = loadDataDir(myLen,Eig0,EIGFILE)
     end if

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveVOSBEig()
     logical :: saveData, saveDataDir

     saveVOSBEig = .false.

     if ( sVOSB <=0 ) return

     if (sVOSB > SEQDIRNUM) then
        saveVOSBEig = saveData(myVOSB%mLen,VOSB,srMode,VOSBFILE)
        if (saveVOSBEig) saveVOSBEig = saveData(myLen,Eig0,srMode,EIGFILE)
     else
        saveVOSBEig = saveDataDir(myVOSB%mLen,VOSB,VOSBFILE)
        if (saveVOSBEig) saveVOSBEig = saveDataDir(myLen,Eig0,EIGFILE)
     end if

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadDEP()
     logical :: loadDataDir, loadDataDir_CX

     loadDEP=.false.
     if (totalDEP) then
        if (sCX) then
           loadDEP = loadDataDir_CX(myDEP%mLen,DEPCX,DEPFILE)
        else
           loadDEP = loadDataDir(myDEP%mLen,DEP,DEPFILE)
        end if
      end if
end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadOutH()
     logical :: loadDataDir, loadDataDir_CX

     loadoutH=.true.
     if (sNDVR) then
        if (sCX) then
           loadOutH = loadDataDir_CX(outLen,OUTHCX,OUTHFILE)
        else
           loadOutH = loadDataDir(outLen,OUTH,OUTHFILE)
        end if
     end if

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadHW()
     integer ::  info, cnt

     loadHW = .false.
     if (srMode) then
          open(99,FILE=HWFile, STATUS='OLD', FORM='UNFORMATTED', &
               IOSTAT=info)     
          if (info /= 0) return
          read(99) cnt
          if (cnt/=sOSBW%mCnt) return
	  read(99,IOSTAT=info) HW                         
          close(99)
      else
          open(99,FILE=HWFile, STATUS='OLD', IOSTAT=info)
          if (info /= 0) return
          read(99,*) cnt
          if (cnt/=sOSBW%mCnt) return
	  read(99,*, IOSTAT=info) HW      
          close(99)
     end if

     loadHW = (info==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadHWCX()
     integer ::  info, cnt

     loadHWCX = .false.
     if (srMode) then
          open(99,FILE=HWFile, STATUS='OLD', FORM='UNFORMATTED', &
               IOSTAT=info)     
          if (info /= 0) return
          read(99) cnt
          if (cnt/=sOSBW%mCnt) return
	  read(99,IOSTAT=info) HWCX                         
          close(99)
      else
          open(99,FILE=HWFile, STATUS='OLD', IOSTAT=info)
          if (info /= 0) return
          read(99,*) cnt
          if (cnt/=sOSBW%mCnt) return
	  read(99,*, IOSTAT=info) HWCX      
          close(99)
     end if

     loadHWCX = (info==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveHW()
    integer :: info
 
    saveHW = .false.

    if (srMode) then
          open(99,FILE=HWFile, STATUS='REPLACE', FORM='UNFORMATTED', &
               IOSTAT=info)
          if (info /= 0) return
          write(99) sOSBW%mCnt
          write(99,IOSTAT=info) HW
          close(99)
      else
          open(99,FILE=HWFile, STATUS='REPLACE', IOSTAT=info)
          if (info /= 0) return
          write(99,*) sOSBW%mCnt
	  write(99,*,IOSTAT=info) HW
          close(99)
     end if
     
     saveHW = (info==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveHWCX()
    integer :: info
 
    saveHWCX = .false.

    if (srMode) then
          open(99,FILE=HWFile, STATUS='REPLACE', FORM='UNFORMATTED', &
               IOSTAT=info)
          if (info /= 0) return
          write(99) sOSBW%mCnt
          write(99,IOSTAT=info) HWCX
          close(99)
      else
          open(99,FILE=HWFile, STATUS='REPLACE', IOSTAT=info)
          if (info /= 0) return
          write(99,*) sOSBW%mCnt
	  write(99,*,IOSTAT=info) HWCX
          close(99)
     end if
     
     saveHWCX = (info==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

