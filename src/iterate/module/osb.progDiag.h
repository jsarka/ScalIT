!cccccccccccccccccccccccccccccccccccccccccccccc
!c      subroutine to do HOSB_DIAG            c
!cccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function progDiag( )

    select case (sVOSB)
    case (0:)    
         call doDiag()
         call saveDiagData()
         progDiag = .true.

    case (:-1)
         progDiag = loadDiagData()

    end select

    if (sOSB==TOSBW) then
        call adjustWinSize(sOSBW%mCnt, sOSBW%mE0, sOSBW%mDE)
    else
        sOSBW%mCnt=getOSBWSize(sOSBW%mE0,sOSBW%mDE)
    end if

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine doDiag( )

    double precision :: ct1, ct2

    print *
    print *, 'ccccccccccccccccccccccccccccccccccccccccccccccccc'
    print *, 'c        Diagonalization of HOSB                c'
    print *, 'ccccccccccccccccccccccccccccccccccccccccccccccccc'
    print *      

    print *, '  Diagonalizing HOSB ......'

    call CPU_TIME(ct1)

    if (.NOT. BlockDiag())    &
       print *, ' Not Efficient Number of Jacobi Iteration in Block Diagonalization'

    call CPU_TIME(ct2)

    print *, 'Time for HOSB Diagonization(sec):',ct2-ct1
    print *, '=============    End of Diagonalization ========='

end  subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine saveDiagData()
   double precision :: time1, time2

   print *
   print *, 'cccccccccccccccccccccccccccccccccccccccccccccccc'
   print *, 'c       Save HOSB, VOSB, EIGVal Results        c'
   print *, 'cccccccccccccccccccccccccccccccccccccccccccccccc' 
   print *
   print *, '  Intemediate Data(HOSB, VOSB, Eig0) are saved in Binary file'

   if (sHOSB > 0 ) then
      print *
      print *, '  Saving HOSB ....... '
      write( *,10) ' HOSB is saved in file: ', HOSBFile
      if (sHOSB>SEQDIRNUM) then
         print *, ' HOSB data are stored in sequential access mode'
      else
         print *, ' HOSB data are stored in direct access mode'
      end if

      if (sST) then     
         call CPU_Time(time1)
         if ( .NOT. saveHOSB() )   print *, ' Error in Saving HOSB!'  
         call CPU_TIME(time2)
         print *, '  Time to save HOSB:(sec)',(time2-time1) 
      else
         print *, '  HOSB are stored during Block Jacobi Diagonalization!'
      end if
   end if

   if (sVOSB > 0) then
      print *
      print *, '  Saving VOSB/Eig ....... '
      write( *,10) ' VOSB is saved in file: ', VOSBFile
      write( *,10) ' Eig0 is saved in file: ', EigFile
      if (sVOSB>SEQDIRNUM) then
         print *, ' VOSB/Eig0 data are stored in sequential access mode'
      else
         print *, ' VOSB/Eig0 data are stored in direct access mode'
      end if

      call CPU_Time(time1)
      if ( .NOT. saveVOSBEig()) print *, ' Error in Savine VOSB/Eig0!'
      call CPU_TIME(time2)
      print *, 'Time to save VOSB/Eig0:(sec)',(time2-time1) 
   end if

   print *
   print *, '========  Finish HOSB/VOSB/EIGVAL  Saving  ======='
   print *
   print *  

   10 FORMAT(2x, A, A)

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadDiagData()
   double precision :: time1, time2

   loadDiagData = .false.

   print *
   print *, 'cccccccccccccccccccccccccccccccccccccccccccccccc'
   print *, 'c         Read   HOSB,   VOSB,   EIGVal       c'
   print *, 'ccccccccccccccccccccccccccccccccccccccccccccccc' 
   print *
   print *, '  Intemediate Data (HOSB, VOSB, Eig) are loaded from Binary File'

   if ((sHOSB<0)) then
      print *
      print *, ' Loading HOSB data ...... '
      if (sST) then     
         write(*,10), ' Loading HOSB from file: ', HOSBFILE
        if (sHOSB<-SEQDIRNUM) then
           print *, ' HOSB data are loaded in sequential access mode'
        else
           print *, ' HOSB data are loaded in direct access mode'
        end if

         call CPU_Time(time1)
         if (.NOT. loadHOSB()) then
            print *, ' Error in loading HOSB !'
            return
         end if
         call CPU_TIME(time2)
         print *, 'Time to load HOSB:(sec)',(time2-time1) 
      else
         print *, ' Partial HOSB are allocated, HOSB is not loaded!'
      end if
   end if

   if (sVOSB<0) then
      print *
      print *, ' Loading VOSB/Eig0 data ....... '
      write(*,10) ' Loading VOSB from file: ', VOSBFILE
      write(*,10) ' Loading Eig0 from file: ', EIGFILE
      if (sVOSB<-SEQDIRNUM) then
         print *, ' VOSB/Eig0 data are loaded in sequential access mode'
      else
         print *, ' VOSB/Eig0 data are loaded in direct access mode'
      end if

      call CPU_Time(time1)
      if ( .NOT. loadVOSBEig()) then
          print *, ' Error in loading VOSB/Eig0 ! '
          return
      end if
      call CPU_TIME(time2)
      print *, 'Time to load VOSB/Eig0:(sec)',(time2-time1) 
   end if
  
   print *, '***********  Finish HOSB/VOSB/EIGVAL  Reading  ************'
   print *

   loadDiagData = .true.
   
   10 FORMAT(2x, A, A)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
