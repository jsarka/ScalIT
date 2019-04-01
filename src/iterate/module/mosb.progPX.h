!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            PX Testing Subroutine                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine  progPX()

    double precision, dimension(plen(sF)) :: X0, X
    double complex,   dimension(plen(sF)) :: X0CX, XCX    
    double precision :: t1, t2
    logical :: log0

    if (id==rootID) then
       print *
       print *, '*****************************************'
       print *, '*        Testing  PX method             *'
       print *, '*****************************************'
       print *       

       if (sCX) then
          print *, '    P*X: P is real, X is complex'
       else
          print *, '    P*X: P is real, X is real'
       end if

       select case (sOSB)
       case (TOSB)
           print *, '    OSB Preconditioner'
       case (TOSBD1)
           print *, '    OSBD1 Preconditioner'
       case (TOSBD2)
           print *, '    OSBD2 Preconditioner'
       case (TOSBW)
           print *, '    OSBW Preconditioner'
       end select

       select case (sPC)
       case (:TA0)
           print *, '    Calculating [1/H]*X'
       case (TA1)
           print *, '    Calculating [1/(H-E)]*X'
       case (TAAP)
           print *, '    Calculating [1/(H-E+iAP)]*X'
       case (TAAPP)
           print *, '    Calculating [1/(H-E+iAPP)]*X'
       case (TAAPR)
           print *, '    Calculating [1/(H-E+iAPR)]*X'
       case (TMAP)
           print *, '    Calculating [1/(H-E-iAP)]*X'
       case (TMAPP)
           print *, '    Calculating [1/(H-E-iAPP)]*X'
       case (TMAPR)
           print *, '    Calculating [1/(H-E-iAPR)]*X'
       end select

       t1 = MPI_WTime()

    end if

    if (sOSB==TOSBW)  log0 = initMOSBW()

    if (sCX)  then   
         call randVec_cx(plen(sF), X0CX)
         call PX_DX(plen(sF),X0CX, XCX)
    else
         if (sAP) then
            call randVec_cx(plen(sF), X0CX)
            call PX_DX(plen(sF),X0CX, XCX)
         else
            call random_number(X0)
            call PX(plen(sF),X0, X)
         end if
    end if

    if (sOSB==TOSBW)  call finalMOSBW()

    if (id==rootID) then
       t2 = MPI_WTime()   

       print *, '  Time to calculate P*X:', t2-t1
       print *, '*****************************************'
       print *
    end if
 end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine  progPXFile(x0fname, svMode, x1fname, pxfname)
    character(len=*), intent(IN) :: x0fname,x1fname, pxfname
    logical, intent(IN) :: svMode

    double precision, dimension(plen(sF)) :: X0, X
    double complex,   dimension(plen(sF)) :: X0CX, XCX    
    double precision :: t1, t2
    logical :: log0

    if (id==rootID) then
       print *
       print *, '*****************************************'
       print *, '*        Testing  PX method             *'
       print *, '*****************************************'
       print *       


       print *, ' Read initial X to file:', x0fname
       if (sCX) then
          print *, ' P*X: P is real, X is complex'
       else
          print *, ' P*X: P is real, X is real'
       end if
       print *, ' Save initial X to file:', x1fname
       print *, ' Save final H*X to file:', pxfname

       select case (sOSB)
       case (TOSB)
           print *, '    OSB Preconditioner'
       case (TOSBD1)
           print *, '    OSBD1 Preconditioner'
       case (TOSBD2)
           print *, '    OSBD2 Preconditioner'
       case (TOSBW)
           print *, '    OSBW Preconditioner'
       end select

       select case (sPC)
       case (:TA0)
           print *, '    Calculating [1/H]*X'
       case (TA1)
           print *, '    Calculating [1/(H-E)]*X'
       case (TAAP)
           print *, '    Calculating [1/(H-E+iAP)]*X'
       case (TAAPP)
           print *, '    Calculating [1/(H-E+iAPP)]*X'
       case (TAAPR)
           print *, '    Calculating [1/(H-E+iAPR)]*X'
       case (TMAP)
           print *, '    Calculating [1/(H-E-iAP)]*X'
       case (TMAPP)
           print *, '    Calculating [1/(H-E-iAPP)]*X'
       case (TMAPR)
           print *, '    Calculating [1/(H-E-iAPR)]*X'
       end select

       t1 = MPI_WTime()

    end if

    if (sOSB==TOSBW)  log0 = initMOSBW()

    if (sCX)  then   
         log0 = loadDiag_CX(.true.,x0fname,X0CX)
         log0 = saveDiag_CX(svMode,x1fname,X0CX)
         call PX_DX(plen(sF),X0CX, XCX)
         log0 = saveDiag_CX(svMode,pxfname,XCX)
    else
         if (sAP) then
            log0 = loadDiag_CX(.true.,x0fname,X0CX)
            log0 = saveDiag_CX(svMode,x1fname,X0CX)
            call PX_DX(plen(sF),X0CX, XCX)
            log0 = saveDiag_CX(svMode,pxfname,XCX)
         else
            log0 = loadDiag(.true.,x0fname,X0)
            log0 = saveDiag(svMode,x1fname,X0)
            call PX(plen(sF),X0, X)
            log0 = saveDiag(svMode,pxfname,X)
         end if
    end if

    if (sOSB==TOSBW)  call finalMOSBW()

    if (id==rootID) then
       t2 = MPI_WTime()

       print *, '  Time to calculate P*X:', t2-t1
       print *, '*****************************************'
    end if
 end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
