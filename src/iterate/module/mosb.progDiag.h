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

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine doDiag( )

    double precision :: ct1, ct2, MPI_WTime

    if (id==rootID) then
        print *, 'ccccccccccccccccccccccccccccccccccccccccccccccccc'
        print *, 'c        Diagonalization of HOSB                c'
        print *, 'ccccccccccccccccccccccccccccccccccccccccccccccccc'
        print *, '  Diagonalizing HOSB ......'
        ct1 = MPI_WTime()
    end if

    if (.NOT. BlockDiag())  then
       if (id==rootID) &
          print *, ' Not Efficient Number of Jacobi Iteration in Block Diagonalization'
    end if

    if (id==rootID) then
        ct2 = MPI_WTIME()
        print *
        print *, '    ---------------------------------------------------------------------'
        print *, '      MPI WTime for HOSB Diagonization(sec):',ct2-ct1
        print *, '    ---------------------------------------------------------------------'
        print *
        print *, '=============    End of Diagonalization ========='
    end if

end  subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine saveDiagData()
   double precision :: time1, time2, MPI_WTime

   if (id==rootID)then
       print *
       print *, 'cccccccccccccccccccccccccccccccccccccccccccccccc'
       print *, 'c       Save HOSB, VOSB, EIGVal Results        c'
       print *, 'cccccccccccccccccccccccccccccccccccccccccccccccc' 
       print *, '  Intemediate Data(HOSB, VOSB, Eig0) are saved in Binary file'
   end if
 
   if ((sHOSB>0) ) then
      if (id==rootID) then
          print *
          print *, '  Saving HOSB ....... '
          write( *,10) ' HOSB is saved in file: ', fHOSB
          print *, ' HOSB data are stored in MPI direct access mode'
          time1 = MPI_WTime()
      end if

      if (sST) then
         if ( .NOT. saveHOSB() ) then
             if (id==rootID)  print *, ' Error in Saving HOSB!'
         end if
      else
         if (id==rootID)  print *, ' HOSB is stored during Block Jacobi!'
      end if

      if (id==rootID) then      
         time2 = MPI_WTIME()
         print *, '  MPI WTime to save HOSB:(sec)',(time2-time1) 
      end if
   else
      if (id==rootID)  print *, ' HOSB is not stored!'
   end if


   if ( sVOSB > 0) then
      if (id==rootID) then
          print *
          print *, '  Saving VOSB/Eig ....... '
          write( *,10) ' VOSB is saved in file: ', fVOSB
          write( *,10) ' Eig0 is saved in file: ', fEig
          print *, ' VOSB/Eig0 data are stored in MPI direct access mode'

          time1 = MPI_WTime()
      end if

      if ( .NOT. saveVOSBEig() ) then
         if (id==rootID)  print *, ' Error in Savine VOSB/Eig0!'
      end if

      if (id==rootID ) then
          time2 =  MPI_WTIME()
          print *, 'Time to save VOSB/Eig0:(sec)',(time2-time1) 
      end if
   else
      if (id==rootID)  print *, ' VOSB/Eig0 are not stored!'
   end if

   if (id==rootID ) then
      print *
      print *, '========  Finish HOSB/VOSB/EIGVAL  Saving  ======='
      print *
   end if

   10 FORMAT(2x, A, A)

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadDiagData()
   double precision :: time1, time2

   loadDiagData = .false.

   if (id==rootID ) then
      print *
      print *, 'cccccccccccccccccccccccccccccccccccccccccccccccc'
      print *, 'c         Read   HOSB,   VOSB,   EIGVal       c'
      print *, 'ccccccccccccccccccccccccccccccccccccccccccccccc' 
      print *
      print *, '  Intemediate Data (HOSB, VOSB, Eig) are loaded from Binary File'
   end if

   if ( sHOSB<0 ) then
      if (id==rootID) then
         print *
         print *, ' Loading HOSB data ...... '
         write(*,10), ' Loading HOSB from file: ', fHOSB
         print *, ' HOSB data are loaded in MPI direct access mode!'
         time1 = MPI_WTime()
      end if

      if (sST) then 
         if (.NOT. loadHOSB()) then
             if (id==rootID)  print *, ' Error in loading HOSB !'          
         end if    
      else
         if (id==rootID)  print *, ' HOSB is not all stored!'
      end if  

      if (id==rootID) then
         time2 = MPI_WTIME()
         print *, 'Time to load HOSB:(sec)',(time2-time1) 
      end if
   else
      if (id==rootID)  print *, ' HOSB is not loaded!'
   end if


   if (sVOSB<0) then
      if (id==rootID) then
          print * 
          print *, ' Loading VOSB/Eig0 data ....... '
          write(*,10) ' Loading VOSB from file: ', fVOSB
          write(*,10) ' Loading Eig0 from file: ', fEIG
          print *, ' VOSB/Eig0 data are loaded in MPI direct access mode'
          time1 = MPI_WTime()
      end if

      loadDiagData = loadVOSBEig()

      if ( .NOT. loadDiagData) then
          if (id==rootID) print *, ' Error in loading VOSB/Eig0 ! '
      end if

      if (id==rootID) then
         time2 = MPI_WTIME()
         print *, 'Time to load VOSB/Eig0:(sec)',(time2-time1) 
      end if
   else
      if (id==rootID)  print *, ' VOSB/Eig0 are not loaded!'
   end if

   if (id==rootID ) then
      print *
      print *, '========  Finish HOSB/VOSB/EIGVAL Loading  ======='
      print *
   end if
   
   10 FORMAT(2x, A, A)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
