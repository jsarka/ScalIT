!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine readOSBSTD()
     call myreadOSB(STDFH)
end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine readOSBFile(filename)
     character(len=*), intent(IN) :: filename
     
     open(99, File=filename, status='OLD')
     call myreadOSB(99)
     close(99)    

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Print information about OSB          c
!cccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine printOSB()

      integer :: I
      logical :: myDebug = .false.   ! True will print data size information

      print *
      print *, 'ccccccccccccccccccccccccccccccccccccccccccc'
      print *, 'c   Information about OSB Configuration   c'
      print *, 'ccccccccccccccccccccccccccccccccccccccccccc'
      print *
      print *, '==========================================='
      print *, '             Basic Information             '
      print *, '==========================================='
      print *, ' Number of Layers:', sF
      print *, ' Number of Points for Each Layer:'
      print *, (sN(I), i=1, sF)
      print *, ' Coordinate Dependency:', totalDep
      print *, ' Coord. Dep. at each Layer:',(sDEP(I),i=1,sF)
      print *

      select case (sJOB)
      case (JOB_RES1,JOB_RES2)
          print *, ' Calculate the Resonance States using PIST algorithm!'
          write(*,190) 'File name for absorption potential:', APPFile
          print *
          print *, ' Parameters for PIST/LANCZOS Convergence'
          write(*,125) sConv%mE0, sConv%mTol 
          write(*,120) sConv%mNum, sConv%mGap
          write(*,110) sConv%mStart, sConv%mStep, sConv%mMax


      case (JOB_CRP, JOB_CRP1, JOB_CRP2)
          print *, ' Calculate the CRP!'
          write(*,190) 'File name for product absorption potential:', APPFile
          write(*,190) 'File name for reactant absorption potential:', APRFile
          write(*,126) sConv%mE0, sConv%mTol
          write(*,127) sConv%mNum,sConv%mGap
      case default
          print *, ' Calculate the Bound States using PIST algorithm!'
          print *, ' No absorption potential for the calculation'
          print *
          print *, ' Parameters for PIST/LANCZOS Convergence'        
          write(*,125) sConv%mE0, sConv%mTol 
          write(*,120) sConv%mNum, sConv%mGap
          write(*,110) sConv%mStart, sConv%mStep, sConv%mMax

      end select

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
           print *, ' Using simple OSB PreConditioner: P=(Eig-E0)^-1'
      case (TOSBD1)
           print *, ' Using the first OSBD Preconditioner:'
           print *, ' P=(Eig-E0)^-1 or P=mDE^-1'
           write(*,105) sOSBW%mE0,  sOSBW%mDE
      case (TOSBD2)
           print *, ' Using the second OSBD Preconditioner:'
           print *, ' P=(Eig-E0)^-1 or P=[alpha+beta*(Eig-E0)]^-1'
           write(*,106) sOSBW%mE0, sOSBW%mDE, sOSBW%mBeta
      case (TOSBW)
           print *, ' Using the full OSBW Preconditioner:'
           print *, ' P=(Eig-E0)^-1 or P=Hp^-1'
           write(*,107) sOSBW%mE0,  sOSBW%mCnt
      end select

      write(*,130) sCX, (.NOT.sNDVR)
      write(*,135) sST,  sAP
      write(*,140) sHij, totalDep

      print *
      print *, ' Control Paramerters to Save/Load HOSB, VOSB, EIG, HW, VX, PT'
      print *, ' 0:No Save/Load            1: Save         -1:Load'
      write(*, 160) sHOSB, sVOSB, sHW, sVX, sPT

      if (mydebug) then
         print *      
         print *, ' ==========================================================='
         print *, '                    Data Structure in Memory       '
         print *, ' ==========================================================='
         print *
         print *, ' Data Structure for HO '
         call  printDataInfo(myH0)

         print *
         print *, ' Data Structure for HOSB '
         call  printDataInfo(myHOSB)

         print *
         print *, ' Data Structure for VOSB'
         call printDataInfo(myVOSB)

         print *
         print *, ' Data Structure for H0DEP '
         call  printDataInfo(myDEP)

         print *
         print *,  ' General Size of data'
         write(*,170) myLen, outLen, hwLen

         print *
         print *, '========================================================='
         print *, ' Layer | Grid Size |   Dimension  |   # Blocks |  Total  '
         print *, '========================================================='
         do i=1, sF
            write(*,180) i, sN(i), myDim(i), myBlk(i), sN(i)*myDim(i)*myBlk(i)
         end do
      end if
 
      print *
      write(*,190) 'File Name for H0:', H0File
      if (sNDVR) then
         write(*,190) 'File Name for OUTH:', OUTHFile
      else
         write(*,190) 'File Name for RES:', RESFile
      end if
      
      do i = 1, SF
         if (sDep(i)) then
             write(*,190) 'File Name for DEP:', DEPFile
             exit
         end if
      end do


      if (sHOSB /= 0) write(*,190) 'File Name for H0SB:', HOSBFile
      if (sVOSB /= 0) then
          write(*,190) 'File Name for VOSB:', VOSBFile
          write(*,190) 'File Name for EIG0:', EIGFile
      end if
      if (sHW   /= 0) write(*,190) 'File Name for HW:', HWFile
      if (sVX   /= 0) write(*,190) 'File Name for VX:', VXFile
      if (sPT   /= 0) write(*,190) 'File Name for PT:', PTFile

      print *, '===========  End of OSB Parameters  =============='
      print *
    
      100 format('  Max. # of Iter.:',I10, '.    Conv. Tol.', F15.9)
      105 format('  Central Energy E0=',F15.9, '.   Window Width mDE=',F15.9)
      106 format('  Central Energy E0=',F15.9, '.   Window Width mDE=',F15.9, &
                 '.   Beta=', F15.9)
      107 format('  Central Energy E0=',F15.9, '.   Window Size mCnt=',I10)
      110 format('  Start Step:',I6, '.  Increment Step:',I6, '.   Maximum Step:',I6)
      120 format('  # of Interested States:',I5, '.  Gap Step to Compared:',I5)
      125 format('  Central Energy:',F15.9, '.   Conv. Tol.:',F15.9)
      126 format('  Starting Energy:',F15.9, '.   Energy step.:',F15.9)
      127 format('  Number of Eigenvalues:',I5, '.   Number of Energy:',I5)
      130 format('  Complex Version:', L5, '.   Out-most layer is DVR: ', L5)
      135 format('  All HOSB is in Memory:', L5,                           &
                 '. Preconditioner Using Absorption Potential:', L5)
      140 format('  CalHij used in OSBW:', I5,  '.   Has Coordinate Dependency:',L5)

      160 format('  sHOSB:',I2,'.   sVOSB/sEig:',I2,'.  sHW:',I2,'.   sVX:',I2,'.  sPT:',I2)
      170 format('  Total Size:',I10,'.   Out-most Layer Size:',I10,'.   Windows Size:',I5)
      180 format(2X, I4, 4x, I4, 4X, 3(I10,4x))
      190 format(3X, A, A)

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine myReadOSB(fd)
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

     read(fd, '(A)') H0File

     if (sNDVR) then
        read(fd, '(A)') OUTHFile
     else
        read(fd, '(A)') RESFile
     end if
     
     do i = 1, sF
        if (sDep(i)) then
            read(fd, '(A)') DEPFile;      exit
        end if
     end do

     select case (sJOB)
     case (JOB_RES1, JOB_RES2)
        read(fd,'(A)') APPFile    
     case (JOB_CRP,JOB_CRP1,JOB_CRP2)
        read(fd, '(A)') APPFile
        read(fd, '(A)') APRFile
     end select

     if (sHOSB /= 0) read(fd, '(A)') HOSBFile     
     if (sVOSB /= 0) then
         read(fd, '(A)') VOSBFile
         read(fd, '(A)') EigFile
     end if
     if (sHW   /= 0) read(fd, '(A)') HWFile
     if (sVX   /= 0) read(fd, '(A)') VXFile
     if (sPT   /= 0) read(fd, '(A)') PTFile

end  subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
