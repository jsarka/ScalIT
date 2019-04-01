!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                HX Testing Subroutine                  c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine progPX()
  
    double precision :: t0, t1, t2, t3
    double precision :: X0(myLen), X(myLen)
    double complex :: Y0(myLen), Y(myLen)

    print *
    print *, '*****************************************'
    print *, '*        Testing  PX method             *'
    print *, '*****************************************'
    print *
    print *, ' P*X mode:', sPC
    print *, ' Preconditioner mode:', sOSB
    print *, ' Vector Length:', myLen

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
           call PX_DX(myLen, Y0, Y)
        else
           if (sJOB==JOB_BOUND) then
               call random_number(X0)
               call PX(myLen,X0,X)
           else
               call randMat_cx(myLen, 1, Y0)
               call PX_DX(myLen, Y0, Y)
           end if
           call CPU_Time(t3)
           print *, '  Time to calculate P*X:', t3-t2           
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



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine progPXFile(saveMode1, xfname, saveMode, pxfname)
     character(len=*),intent(IN) :: xfname,pxfname
    logical, intent(IN) :: saveMode1, saveMode 
  
    double precision :: t0, t1, t2
    double precision :: X0(myLen), X(myLen)
    double complex :: Y0(myLen), Y(myLen)
    logical :: saveData, saveData_CX  

    print *
    print *, '*****************************************'
    print *, '*        Testing  PX method             *'
    print *, '*****************************************'
    print *
    print *, ' P*X mode:', sPC
    print *, ' Preconditioner mode:', sOSB
    print *, ' Vector Length:', myLen,                 &
             ' Save Mode:(T=binary, F=ASCII):', saveMode
    print *, ' The initial vector is stored in file:', xfname
    print *, ' The final P*X vector is stored in file:', pxfname

    call CPU_TIME(t0)

    print *
    if (sOSB==TOSBW) then
       if (.NOT.initOSBW()) then 
           print *, ' Error in initalizing OSBW '   
           return
       end if
       call CPU_Time(t1)
       print *, ' Time to calculate Hij for Energy Window:',t1-t0
    end if

    if (sCX)  then   
       call randMat_cx(myLen, 1, Y0)
       call PX_DX(myLen, Y0, Y)

       if (.NOT. saveData_CX(myLen,Y0,saveMode1, xfname))   &
           print *, ' Error in saving initial vector in file.'

       if (.NOT. saveData_CX(myLen,Y,saveMode, pxfname))   &
          print *, ' Error in saving final vector P*X in file.'
    else
       if (sJOB==JOB_BOUND) then
           call random_number(X0)
           call PX(myLen,X0,X)

           if (.NOT. saveData(myLen,X0,saveMode1, xfname))   &
               print *, ' Error in saving initial vector in file.'
           if (.NOT. saveData(myLen,X,saveMode, pxfname))    &
               print *, ' Error in saving final vector P*X in file.'
          !print *, 'initial vector:', X0
          !print *, 'P*X vector:',X
       else
           call randMat_cx(myLen, 1, Y0)
           call PX_DX(myLen, Y0, Y)

           if (.NOT. saveData_CX(myLen,Y0,saveMode1, xfname))   &
               print *, ' Error in saving initial vector in file.'
           if (.NOT. saveData_CX(myLen,Y, saveMode, pxfname))   &
               print *, ' Error in saving final vector P*X in file.'
       end if
     
       call CPU_Time(t2)
       print *, '  Time to calculate P*X:', t2-t1           
    end if
 
    print *, '*****************************************'
    print *

 end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

