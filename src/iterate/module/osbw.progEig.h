!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         Eigen Value Testing Subroutine              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine progPistEig()

    double precision :: eig(sConv%mNum),crp(sConv%mNum) 
    double complex   :: eigCX(sConv%mNum)

    integer :: iter, i
    double precision  :: ct1, ct2, maxRES

    call CPU_TIME(ct1)
    print *
    print *, '================================================'
    print *, '     Get  Eigen Values Using PIST method        '
    print *, '================================================'

    call printOSBWInfo()

    if (sVX<0)  print *, '  Final Vector are stored in file:', VXFILE

    select case (sJOB)
    case (JOB_RES1)    ! resounace states
        if (sCX) then
	        print *, ' ------ Resonance states for Complex Version of H0  -------'
	        print *, ' ------ Using Symmetric version of PIST  -------'
           iter = OSB_PISTCONVCX1(sConv%mNum, EigCX, maxRes)        
           print *
           print *,'Max Lanczos Error:', maxres, ' Eigen Values:'
           print 1001, eigCX(1:sConv%mNum)
        else
           print *, ' ------  Resonance states for Real Version of H0  -------'
	   print *, ' ------ Using Symmetric version of PIST  -------'
           iter = OSB_PISTCONVDX1(sConv%mNum, EigCX, maxRes)
           print * 
           print *,'Max Lanczos Error:', maxres, ' Eigen Values:'
           print 1001, eigCX(1:sConv%mNum) 
       end if

    case (JOB_RES2)    ! resounace states
        if (sCX) then
	        print *, ' ------ Resonance states for Complex Version of H0 -------'
	        print *, ' ------ Using Conjugate version of PIST  -------'
           iter = OSB_PISTCONVCX2(sConv%mNum, EigCX, maxRes)        
           print *
           print *,'Max Lanczos Error:', maxres, ' Eigen Values:'
           print 1001, eigCX(1:sConv%mNum)
        else
           print *, ' ------  Resonance states for Real Version of H0:Version 2  -------'
	        print *, ' ------ Using Conjugate version of PIST  -------'
           iter = OSB_PISTCONVDX2(sConv%mNum, EigCX, maxRes)
           print * 
           print *,'Max Lanczos Error:', maxres, ' Eigen Values:'
           print 1001, eigCX(1:sConv%mNum) 
       end if

    case (JOB_CRP,JOB_CRP1,JOB_CRP2)    ! CRP
        if (sCX) then
	        print *, ' ------ CRP calculation for Complex Version of H0  -------'
        else
           print *, ' ------  CRP calculation for Real Version of H0  -------'
        end if

   	do i = 1, sConv%mNum
  	   eig(i)=sConv%mE0 + (i-1)*sConv%mTol
        end do

	call osbCRP(sConv%mNum, sConv%mGap, eig, crp)

    case default   ! bound states
        if (sCX) then
	        print *, ' ------ Bound states for Complex Version of H0  -------'
           iter = OSB_PISTCONVCX0(sConv%mNum, EigCX, maxRes)        
           print *
           print *,'Max Lanczos Error:', maxres, ' Eigen Values:'
           print 1001, eigCX(1:sConv%mNum)
        else
           print *, ' ------  Bound states for Real Version of H0  -------'
           iter = OSB_PISTCONV(sConv%mNum, Eig, maxRes)    
           print *    
           print *,'Max Lanczos Error:', maxres, ' Eigen Values:'
           print 1000, eig(1:sConv%mNum)
        end if
    end select

    call CPU_TIME(ct2)
    print *
    print *, 'Time to get PIST Eigen values:(sec)',(ct2-ct1)
    print *, '===========    End of PIST Method   ==================='
    print *
    
 1000 format (E25.15,2x,E25.15,2x,E25.15)
 1001 format (E25.15,2x,E25.15)
 1002 format (E16.10E2)
 
end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


subroutine printOSBWInfo()

    print *
    print *, ' ===== Preconditioner Information ===='
    print *, ' Central Energy:', sOSBW%mE0
    print *
    select case (sOSB)  
    case (TOSB)
        print *, ' Using OSB preconditioner.'
        print *, ' Central Energy:', sOSBW%mE0
    case (TOSBD1)
        print *, ' Using OSBD1 preconditioner '
        print *, ' Central Energy:', sOSBW%mE0, '. Window Width:',sOSBW%mDE
        print *, ' Update Window Size. Window Size:',sOSBW%mCnt
    case (TOSBD2)
        print *, ' Using OSBD2 preconditioner.'
        print *, ' Center Energy:', sOSBW%mE0, '. Window Width:',sOSBW%mDE
        print *, ' Update Window Size. Window Size:',sOSBW%mCnt
        print *, ' Beta Parameter(<-1.0):',sOSBW%mBeta
    case (TOSBW)
        print *, ' Using OSBW full preconditioner.'
        print *, ' Central Energy:', sOSBW%mE0, '. Window Size:',sOSBW%mCnt
        print *, ' Update Window Width. Window Width:',sOSBW%mDe
    case default
        print *, ' Using OSB preconditioner.'
        print *, ' Central Energy:', sOSBW%mE0
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

