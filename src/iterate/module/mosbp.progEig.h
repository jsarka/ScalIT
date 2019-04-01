!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         Eigen Value Testing Subroutine              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine progPistEig()

    double precision :: eig(sConv%mNum), crp(sConv%mNum)
    double complex   :: eigCX(sConv%mNum)

    integer :: iter, ierr, i
    double precision  :: ct1, ct2, maxRES

    if (ID==rootID) then
       ct1 = MPI_WTIME()
       print *
       print *, '================================================'
       print *, '     Get  Eigen Values Using PIST method        '
       print *, '================================================'

       call printOSBWInfo()

       if (sVX>0)  print *, '  Final Vector are stored in file:', fVX
    end if
    
    select case (sJOB)
    case (JOB_RES1)    ! resounace states
        if (ID==rootID) then
           if (sCX)  then
	      print *, ' ------ Resonance states for Complex Version of H0  -------'
           else
              print *, ' ------  Resonance states for Real Version of H0  -------'
           end if
       end if

       if (sCX) then
            iter = OSB_PISTCONVCX1(sConv%mNum, EigCX, maxRes)
       else
            iter = OSB_PISTCONVDX1(sConv%mNum, EigCX, maxRes)
       end if

       if (ID==rootID) then
           print * 
           print *,'Max Lanczos Error', maxres, ' Eigen Values:'
           print 1001, eigCX(1:sConv%mNum) 
       end if

    case (JOB_RES2)    ! resounace states
        if (ID==rootID) then
           if (sCX)  then
	      print *, ' ------ Resonance states for Complex Version of H0  -------'
           else
              print *, ' ------  Resonance states for Real Version of H0  -------'
           end if
       end if

       if (sCX) then
            iter = OSB_PISTCONVCX2(sConv%mNum, EigCX, maxRes)
       else
            iter = OSB_PISTCONVDX2(sConv%mNum, EigCX, maxRes)
       end if

       if (ID==rootID) then
           print * 
           print *,'Max Lanczos Error', maxres, ' Eigen Values:'
           print 1001, eigCX(1:sConv%mNum) 
       end if

    case (JOB_CRP,JOB_CRP1,JOB_CRP2)    ! CRP
        if (ID==rootID) then
           if (sCX) then
	      print *, ' ------ CRP calculation for Complex Version of H0  -------'
           else
              print *, ' ------  CRP calculation for Real Version of H0  -------'
           end if
        end if

   	do i = 1, sConv%mNum
  	   eig(i)=sConv%mE0 + (i-1)*sConv%mTol
        end do

	call osbCRP(sConv%mGap, sConv%mNum, eig, crp)

       if (ID==rootID) then
           print *
           print *,' Energy:                  CRP:'
           do i = 1, sConv%mNum
                print *, eig(i), crp(i) 
           end do
       end if

    case default   ! bound states
       if (ID==rootID) then
           if (sCX)  then
              print *, ' ------ Bound states for Complex Version of H0  -------'
           else
              print *, ' ------  Bound states for Real Version of H0  -------'
              call printConvInfo()
           end if
       end if

       if (sCX) then
            iter = OSB_PISTCONVCX0(sConv%mNum, EigCX, maxRes)
       else
            iter = OSB_PISTCONV(sConv%mNum, Eig, maxRes)
       end if

       if (ID==rootID) then
          if (sCX) then
               print * 
               print *,'Max Lanczos Error', maxres, ' Eigen Values:'
               print 1001, eigcx(1:sConv%mNum)
          else
               print * 
               print *,'Max Lanczos Error', maxres, ' Eigen Values:'
               print 1000, eig(1:sConv%mNum)
               print *
               print *,' eigenvalues:'
               print 1002, eig(1:sConv%mNum)
          end if
       end if

    end select

    if (ID==rootID) then

        ct2 = MPI_WTIME()
        print *
        print *, '   ---------------------------------------------------------------------'
        print *, '     Timing for PIST and Lanczos Diagonalization (MPI_WTime):'
        print *, '        ',(ct2-ct1),'seconds'
        print *, '   ---------------------------------------------------------------------'
        print * 
        print *, '=======================    End of program   ========================='
        print *

    end if

 1000 format (E25.15,2x,E25.15,2x,E25.15)
 1001 format (E25.15,2x,E25.15)
 1002 format (F15.10)

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine printOSBWInfo()

    print *
    print *, ' ===== Preconditioner Information ===='
    print *, ' Central Energy:', sOSBW%mE0
    print *
    select case (sOSB)  
    case (TOSB)
        print *, ' Using OSB preconditioner.'
    case (TOSBD1)
        print *, ' Using OSBD1 preconditioner '
        print *, ' OSB Windows:', sOSBW%mDE
    case (TOSBD2)
        print *, ' Using OSBD2 preconditioner.'
        print *, ' OSB Windows:', sOSBW%mDE,' Beta=',sOSBW%mBeta
    case (TOSBW)
        print *, ' Using OSBW1 full preconditioner.'
        print *, ' OSB Original Window Width:', sOSBW%mDE
        print *, ' # States in window:', sOSBW%mCnt
    case default
        print *, ' Using OSB preconditioner.'
    end select

    print *
    print *, ' ===== End of  Preconditioner Information ===='
    print *
  
end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine printConvInfo()
    print *
    print *, ' ===== Convergence Information ===='
    print *
    print *, ' Central Energy:', sConv%mE0, ' Conv. Tol.:', sConv%mTol
    print *, ' Number of interested states:',sConv%mNum
    print *, ' Start Size:', sConv%mStart, ' Increment Size:', sConv%mStep
    print *, ' Max. Size:', sConv%mMax, ' Conv. Testing Gap:', sConv%mGap
    print *
    print *, ' ===== End of Convergence Information ===='
    print *
  
end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

