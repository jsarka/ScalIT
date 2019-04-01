!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine readMOSBSTD()

     call myreadMOSB(STDFH)

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine readMOSBFile(filename)
     character(len=*), intent(IN) :: filename
     
     open(99, File=filename, status='OLD')
     call myreadMOSB(99)
     close(99)    

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Print information about OSB          c
!cccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine printMOSB()
      integer :: I

      print *
      print *, 'cccccccccccccccccccccccccccccccccccccccccccc'
      print *, 'c   Information about MOSB Configuration   c'
      print *, 'cccccccccccccccccccccccccccccccccccccccccccc'
      print *
      print *, '==========================================='
      print *, '             Basic Information             '
      print *, '==========================================='

      select case (sJOB)
      case (JOB_RES1)
           print *, ' Calculate the Resonance States!'
           print *, ' Using complex Symmetric version of PIST!'
      case (JOB_RES2)
           print *, ' Calculate the Resonance States!'
           print *, ' Using complex Conjugate version of PIST!'

      case (JOB_CRP)
           print *, ' Calculate the CRP!'
           print *, ' Version 0'

      case (JOB_CRP1)
           print *, ' Calculate the CRP!'
           print *, ' Version 1'

      case (JOB_CRP2)
           print *, ' Calculate the CRP!'
           print *, ' Version 2'

      case default
           print *, ' Calculate the Bound States!'
      end select

      print *
      print *, ' Number of Computing Nodes:', nNodes 
      print *, ' ID of current Computing Nodes:', id      
      print *, ' Number of Layers:', sF
      print *, ' Number of Points for Each Layer:'
      print *,  (sN(I), i=1, sF)
      print *, ' Coordinate Dependency:', totalDep
      print *, ' Coord. Dep. at each Layer:',(sDEP(I),i=1,sF)
      print *, ' Splitting Level:', myNode%spLevel
      print *, ' Size for double precision data:', dbSize
      print *, ' Size for double complex data:', cxSize

      print *
      print *, ' Parameters for Block Jacobi Diagaonalization'
      write(*,100)  sBJ%mMax, sBJ%mTol
      print *
      print *, ' Parameters for QMR Linear Solver'
      write(*,100) sQMR%mMax, sQMR%mTol

      print *
      print *, ' Parameters for OSB Preconditioner'
      select case (sOSB)
      case (TOSB)
           print *, ' Simple OSB preconditioner. E0=', sOSBW%mE0
      case (TOSBD1)
           print *, ' First type OSBD preconditioner:'
           write(*,105) sOSBW%mE0,  sOSBW%mDE
      case (TOSBD2)
           print *, ' Second type OSBD preconditioner:'
           write(*,106) sOSBW%mE0,  sOSBW%mDE,sOSBW%mBeta
      case (TOSBW)
           print *, ' OSBW preconditioner:'
           write(*,107) sOSBW%mE0, sOSBW%mCnt
      end select

      print *
      select case (sJOB)
      case (JOB_RES1, JOB_RES2)
           print *, ' Parameters for PIST/LANCZOS Convergence'
           write(*,125) sConv%mE0, sConv%mTol
           write(*,120) sConv%mNum, sConv%mGap
           write(*,110) sConv%mStart, sConv%mStep, sConv%mMax
         
      case (JOB_CRP,JOB_CRP1,JOB_CRP2)
           print *, ' Parameters for CRP Calculation'
           write(*,126) sConv%mE0, sConv%mTol
           write(*,111) sConv%mNum, sConv%mGap
          
      case default
           print *, ' Parameters for PIST/LANCZOS Convergence'
           write(*,125) sConv%mE0, sConv%mTol
           write(*,120) sConv%mNum, sConv%mGap
           write(*,110) sConv%mStart, sConv%mStep, sConv%mMax
           
      end select

      print *
      write(*,130) sCX, .NOT.sNDVR, sST
      write(*,140) sAP, totalDep

      print *
      print *, ' Control Paramerters to Save/Load HOSB, VOSB, EIG, HW, VX, PT'
      print *, '0:No Save/Load            1: Save         -1:Load'
      write(*, 160) sHOSB, sVOSB, sHW, sVX, sPT

      if (mydebug) then
          print *      
          print *, '==========================================================='
          print *, '                    Data Structure in Memory       '
          print *, '==========================================================='

          print *
          print *, ' Global Configure of Data '
          call  printGData(myconf)

          print *
          print *, ' Struture of Current Computing Node '
          call  printMNode(sF, myNode)

          print *
          print *, ' Struture of Data Distribution '
          call  printMData(sF, myData)

          print *
          print *, ' Data Structure for Vi '
          call  printSeqInfo(myVi)

          print *
          print *, ' Data Structure for HO '
          call  printSeqInfo(myH0)

          print *
          print *, ' Data Structure for HOSB '
          call  printMGrid(sF,myHOSB)

          print *
          print *, ' Data Structure for VOSB'
          call printMGrid(sF,myVOSB)

          print *
          print *, ' Data Structure for RES '
          call  printMGrid(sF,myRES)

          print *
          print *, ' Data Structure for DEP '
          call  printMGrid(sF,myDEP)
      end if

      print *
      print *, '===========    Data and Related File name   ========='
      write(*,190) 'File Name for H0:', fH0
      if (sNDVR) then
         write(*,190) 'File Name for OUTH:', fOUTH
      else
         write(*,190) 'File Name for RES:', fRES
      end if
      
      do i = 1, SF
         if (sDep(i)) then
             write(*,190) 'File Name for DEP:', fDEP
             exit
         end if
      end do

      select case (sJOB)
      case (JOB_RES1, JOB_RES2)
         write(*,190) ' File name for absorption potential:', fAPP

      case (JOB_CRP,JOB_CRP1,JOB_CRP2)
         write(*,190) 'File name for product absorption potential:', fAPP
         write(*,190) 'File name for reactant absoprtion potential:', fAPR
      end select 

      if (sHOSB /= 0) write(*,190) 'File Name for H0SB:', fHOSB
      if (sVOSB /= 0) then
          write(*,190) 'File Name for VOSB:', fVOSB
          write(*,190) 'File Name for EIG0:', fEIG
      end if
      if (sHW   /= 0) write(*,190) 'File Name for HW:', fHW
      if (sVX   /= 0) write(*,190) 'File Name for VX:', fVX
      if (sPT   /= 0) write(*,190) 'File Name for PT:', fPT

      print *
      print *, '===========  End of OSB Parameters  =============='
      print *
    
      100 format("  Max. # of Iter.:",I10,  ".    Conv. Tol.", F15.9)
      105 format("  Central Energy:",F15.9, ".  Windows Width:",F15.9)
      106 format("  Central Energy:",F15.9, ".  Windows Width:",F15.9, ".  Beta:", F15.9)
      107 format("  Central Energy:",F15.9, ".  Windows Size:",I10)
      110 format("  Start Step:",I6, '.  Increment Step:',I6, '.   Maximum Step:',I6)
      111 format("  Number of Energy:", I5, ". Number of Eigenvalues:",I5)
      120 format("  # of Interested States:",I5, ".  Gap Step to Compared:",I5)
      125 format("  Central Energy:",F15.9, ".   Conv. Tol.:",F15.9)
      126 format("  Starting Energy:",F15.9, ".  Energy Increment.:",F15.9)
      130 format("  Complex Version:", L5, ".   Out-most layer is DVR: ",L5,        &
                  ".   HOSB is in Memory:", L5)
      140 format("  Use Absorption Potential:",L5, ".   Has Coordinate Dependency:",L5)
      150 format("  H*X choice:",I5, ".   P*X choice:", I5)
      155 format("  Preconditioner:", I5, ". Hij for OSBW: ", I5)
      160 format("  sHOSB:",I2,".   sVOSB:",I2, ".  sHW",I2, ".   sVX:",I2, ".  sPT", I2)
      170 format("  Global Total Size:",I15, ".   Local Total Size:", &
                  I10, ".   Windows Size",I5)
      175 format("  Out-most Layer: Global Total Size:",I15,          &
                 ".   Local Total Size:",I10, ".   Start Position",I15)
      180 format(2X, I4, 2x, 2(I15, 2X), 2(I10,2x))
      185 format(2x, I4, 2X, 3(I15, 2X) )
      190 format(4X, A, A)

end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine myReadMOSB(fd)
     integer, intent(IN) :: fd

     integer :: i

     read(fd, *) sF, (sN(i), i=1, sF)
     read(fd, *) (sDEP(i),i=1,sF);   
     sDEP(1)=.false.

     read(fd, *) sJOB, sOSB

     read(fd, *) sCX, sNDVR, sST, sAP
     if (sNDVR)  sDEP(sF)=.false.

     read(fd, *) sBJ, sQMR

     read(fd, *) sConv

     read(fd, *) sOSBW
     sOSBW%mDE=ABS(sOSBW%mDE)

     read(fd, *) sHOSB, sVOSB, sHW, sVX, sPT     

     read(fd, '(A)') fH0

     if (sNDVR) then
        read(fd, '(A)') fOUTH
     else
        read(fd, '(A)') fRES
     end if
     
     do i = 1, sF
        if (sDep(i)) then
            read(fd, '(A)') fDEP;      exit
        end if
     end do

     select case (sJOB)
     case (JOB_RES1,JOB_RES2)
        read(fd,'(A)') fAPP
     case (JOB_CRP,JOB_CRP1,JOB_CRP2)
        read(fd, '(A)') fAPP
        read(fd, '(A)') fAPR
     end select

     if (sHOSB /= 0) read(fd, '(A)') fHOSB
     if (sVOSB /= 0) then
        read(fd, '(A)') fVOSB
        read(fd, '(A)') fEig
     end if
     if (sHW   /= 0) read(fd, '(A)') fHW
     if (sVX   /= 0) read(fd, '(A)') fVX
     if (sPT   /= 0) read(fd, '(A)') fPT

end  subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine distributeMosb(ierr)
    integer, intent(OUT) :: ierr

    integer, parameter :: NUM_INT=16
    integer, parameter :: NUM_DB=7
    integer, parameter :: NUM_LG=5
    integer :: intData(NUM_INT), i
    logical :: lgData(NUM_LG)
    double precision  :: dbData(NUM_DB)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    if (id==rootID) then
       intData(1)=sF;           
       intData(2)=sBJ%mMax;     intData(3)=sQMR%mMax
       intData(4)=sConv%mStart; intData(5)=sConv%mStep;
       intData(6)=sConv%mMax;   intData(7)=sConv%mNum
       intData(8)=sConv%mGap
       intData(9)=sJOB;         intData(10)=sOSB;
       intData(11)=sHOSB;       intData(12)=sVOSB; 
       intData(13)=sHW;         intData(14)=sVX; 
       intData(15)=sPT;         intData(16)=sOSBW%mCnt; 
    end if
     
    call MPI_BCAST(intData,NUM_INT,MPI_INTEGER,rootID, MPI_COMM_WORLD,ierr)

    sF=intData(1);  
    sBJ%mMax = intData(2);     sQMR%mMax = intData(3)
    sConv%mStart = intData(4); sConv%mStep = intData(5);
    sConv%mMax = intData(6);   sConv%mNum  = intData(7)
    sConv%mGap = intData(8)
    sJOB  = intData(9);   sOSB  = intData(10);
    sHOSB = intData(11);  sVOSB = intData(12)
    sHW   = intData(13);  sVX   = intData(14)
    sPT   = intData(15);  sOSBW%mCnt = intData(16)

    call MPI_BCAST(sN,sF,MPI_INTEGER,rootID, MPI_COMM_WORLD,ierr)

    if (id==rootID) then
       dbData(1)=sBJ%mTol;     dbData(2)=sQMR%mTol
       dbData(3)=sConv%mE0;    dbData(4)=sConv%mTol
       dbData(5)=sOsbw%mE0;    dbData(6)=sOsbw%mDE
       dbData(7)=sOsbw%mBeta
    end if

    call MPI_BCAST(dbData,NUM_DB,MPI_DOUBLE_PRECISION,rootID, MPI_COMM_WORLD,ierr)

    sBJ%mTol  = dbData(1);     sQMR%mTol  = dbData(2)
    sConv%mE0 = dbData(3);     sConv%mTol = dbData(4)
    sOsbw%mE0 = dbData(5);     sOsbw%mDE  = dbData(6)
    sOsbw%mBeta = dbData(7)

    if (id==rootID) then
       lgData(1)=sCX; lgData(2)=sNDVR;  lgData(3)=sST
       lgData(4)=sAP; lgData(5)=srMode
    end if

    call MPI_BCAST(lgData,NUM_LG,MPI_LOGICAL,rootID, MPI_COMM_WORLD,ierr)

    sCX = lgData(1);  sNDVR  = lgData(2);  sST = lgData(3)
    sAP = lgData(4);  srMode = lgData(5)

    call MPI_BCAST(sDep,sF,MPI_LOGICAL,rootID, MPI_COMM_WORLD,ierr)

     totalDep=.false.
     do i = 1, sF
         if (sDep(i)) then
            totalDep=.true.;  exit
         end if
     end do 

    call MPI_BCAST(fH0,MAX_FNAME,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)

    if (sNDVR) then
       call MPI_BCAST(fOUTH,MAX_FNAME,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)
    else
       call MPI_BCAST(fRES,MAX_FNAME,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)
    end if

    if (totalDep)    &
       call MPI_BCAST(fDep,MAX_FNAME,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)
  
    select case (sJOB)
    case (JOB_RES1, JOB_RES2)
       call MPI_BCAST(fAPP,MAX_FNAME,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)

    case (JOB_CRP,JOB_CRP1,JOB_CRP2)
       call MPI_BCAST(fAPP,MAX_FNAME,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(fAPR,MAX_FNAME,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)
    end select

    if (sHOSB/=0)   &
       call MPI_BCAST(fHOSB,MAX_FNAME,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)

    if (sVOSB/=0)  then
       call MPI_BCAST(fVOSB,MAX_FNAME,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(fEIG,MAX_FNAME,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)
    end if

    if (sHW/=0)     &
       call MPI_BCAST(fHW,MAX_FNAME,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)

    if (sVX/=0)     &
       call MPI_BCAST(fVX,MAX_FNAME,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)

    if (sPT/=0)     &
       call MPI_BCAST(fPT,MAX_FNAME,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)

end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
