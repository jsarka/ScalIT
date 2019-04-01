!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                HX Testing Subroutine                  c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine progHX()  
    double precision :: t0, t1, t2
    double precision :: X0(myLen), X(myLen)
    double complex :: Y0(myLen), Y(myLen)

    print *
    print *, '*****************************************'
    print *, '*        Testing  HX method             *'
    print *, '*****************************************'
    print *
    print *, ' H*X mode:', sHC
    print *, ' Vector Length:', myLen

    call CPU_TIME(t0)

    print *, ' Performing Block Jacobi Diagonalization ......'
     if (progDiag()) then
        call CPU_Time(t1)
        print *, '   Time to perform Block Jacobi:', t1-t0
        print *

        if (sCX)  then   
           call randMat_cx(myLen, 1, Y0)
           call HX_CX(myLen, Y0, Y)
        else
           if (sJOB==JOB_BOUND) then
               call random_number(X0)
               call HX(myLen,X0,X)
           else
               call randMat_cx(myLen, 1, Y0)
               call HX_DX(myLen, Y0, Y)
           end if

           print *, '  Time to calculate H*X:', t2-t1           
        end if
    else
        print *, ' Error in Block Jacobi Diagonalization!'
    end if

    call CPU_TIME(t2)   
    print *
    print *, '  Total Time for the calculation:', t2-t0
    print *, '*****************************************'
    print *

 end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine progHXFile(saveMode1, xfname, saveMode, hxfname)
     character(len=*),intent(IN) :: xfname,hxfname
    logical, intent(IN) :: saveMode1, saveMode 

    integer :: i  
    double precision :: t0, t1
    double precision :: X0(myLen), X(myLen)
    double complex :: Y0(myLen), Y(myLen)
    logical :: saveData, saveData_CX  

    print *
    print *, '*****************************************'
    print *, '*        Testing  HX method             *'
    print *, '*****************************************'
    print *
    print *, ' H*X mode:', sHC
    print *, ' Vector Length:', myLen,                 &
             ' Save Mode:(T=binary, F=ASCII):', saveMode
    print *, ' The initial vector is stored in file:', xfname
    print *, ' The final H*X vector is stored in file:',hxfname

    call CPU_TIME(t0)

    print *

    if (sCX)  then   
       call randMat_cx(myLen, 1, Y0)
       call HX_CX(myLen, Y0, Y)

       if (.NOT. saveData_CX(myLen,Y0,saveMode1, xfname))   &
           print *, ' Error in saving initial vectors in file.'
       if (.NOT. saveData_CX(myLen,Y,saveMode, hxfname))   &
           print *, ' Error in saving H*X vectors in file.'
    else
       if (sJOB==JOB_BOUND) then
          call random_number(X0)
          call HX(myLen,X0,X)
          if (.NOT. saveData(myLen,X0,saveMode1, xfname))   &
             print *, ' Error in saving initial vectors in file.'
          if (.NOT. saveData(myLen,X,saveMode, hxfname))   &
             print *, ' Error in saving H*X vectors in file.'
          print *, 'initial vector:', X0
          print *, 'H*X vector:',X 
      else
          call randMat_cx(myLen, 1, Y0)
          call HX_DX(myLen, Y0, Y)

          if (.NOT. saveData_CX(myLen,Y0,saveMode1, xfname))   &
              print *, ' Error in saving initial vectors in file.'
          if (.NOT. saveData_CX(myLen,Y,saveMode, hxfname))   &
              print *, ' Error in saving H*X vectors in file.'
       end if
    end if
 
    call CPU_Time(t1)
         
    print *
    print *, '  Total Time for the calculation:', t1-t0
    print *, '*****************************************'
    print *

 end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


