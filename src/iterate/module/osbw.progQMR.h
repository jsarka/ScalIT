!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                HX Testing Subroutine                  c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine progQMR(fname, saveMode)
    character(len=*),intent(IN) :: fname
    logical, intent(IN) :: saveMode 

    integer :: ncnt  
    double precision :: t0, t1, t2, t3, dre
    double precision :: X0(myLen), X(myLen)
    double complex :: Y0(myLen), Y(myLen)
    logical :: saveData2, saveData2_CX  

    print *
    print *, '*****************************************'
    print *, '*        Testing  HX method             *'
    print *, '*****************************************'
    print *
    print *, ' H*X mode', sPC, ' P*X mode:', sPC
    print *, ' Preconditioner mode:', sOSB
    print *, ' Vector Length:', myLen,                 &
             ' Save Mode:(T=binary, F=ASCII):', saveMode
    print *, ' The initial and final vectors are stored in file:', fname

    call CPU_TIME(t0)

    print *, ' Performing Block Jacobi Diagonalization ......'
     if (progDiag()) then
        call CPU_Time(t1)
        print *, '   Time to perform Block Jacobi:', t1-t0
        print *

        if (sOSB==TOSBW) then
  	   if (.NOT.initOSBW()) then 
 	       print *, ' Error in initalizing OSBW '   
               return
           end if
  	   call CPU_Time(t2)
           print *, ' Time to calculate Hij for Energy Window:',t2-t1
        end if

        if (sCX)  then   
           call randMat_cx(myLen, 1, Y0)
           ncnt = OSB_QMRCX(myLen, Y0, Y, dre)
           if (.NOT. saveData2_CX(myLen,Y0, Y, saveMode, fname))   &
               print *, ' Error in saving vectors in file.'
        else
           if (sJOB==JOB_BOUND) then
               call random_number(X0)
               ncnt = OSB_QMR(myLen, X0, X, dre)
               if (.NOT. saveData2(myLen,X0, X, saveMode, fname))   &
                  print *, ' Error in saving vectors in file.'
           else
               call randMat_cx(myLen, 1, Y0)
               ncnt = OSB_QMRDX(myLen, Y0, Y, dre)
               if (.NOT. saveData2_CX(myLen,Y0, Y, saveMode, fname))   &
                   print *, ' Error in saving vectors in file.'
           end if
           call CPU_Time(t3)
           print *, '  Time to calculate QMR:', t3-t2           
        end if
    else
        print *, ' Error in Block Jacobi Diagonalization!'
    end if

    call CPU_TIME(t3)   
    print *
    print *, '  Total Time for the calculation:', t3-t0
    print *, '*****************************************'
    print *

 end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


