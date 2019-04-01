!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            HX Testing Subroutine                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine  progHX()
    double precision, dimension(plen(sF)) :: X0, X
    double complex,   dimension(plen(sF)) :: X0CX, XCX    
    double precision :: t1, t2

    if (id==rootID) then
       print *
       print *, '*****************************************'
       print *, '*        Testing  HX method             *'
       print *, '*****************************************'
       print *       

       if (sCX) then
          print *, ' H*X: H is complex, X is complex'
       else
          if (sAP) then
             print *, '    H*X: H is real, X is complex'
          else
             print *, '    H*X: H is real, X is real'
          end if
       end if

       select case (sHC)
       case (:TA0)
           print *, '    Calculating H*X'
       case (TA1)
           print *, '    Calculating (H-E)*X'
       case (TAAP)
           print *, '    Calculating (H-E+iAP)*X'
       case (TAAPP)
           print *, '    Calculating (H-E+iAPP)*X'
       case (TAAPR)
           print *, '    Calculating (H-E+iAPR)*X'
       case (TMAP)
           print *, '    Calculating (H-E-iAP)*X'
       case (TMAPP)
           print *, '    Calculating (H-E-iAPP)*X'
       case (TMAPR)
           print *, '    Calculating (H-E-iAPR)*X'
       end select

       t1=MPI_WTime()

    end if

    if (sCX)  then   
         call randVec_cx(plen(sF), X0CX)
         call HX_CX(plen(sF),X0CX, XCX)
    else
         if (sAP) then
             call randVec_cx(plen(sF), X0CX)  
             call HX_DX(plen(sF),X0CX, XCX)
         else
             call random_number(X0)
             call HX(plen(sF),X0, X)
         end if
    end if

    if (id==rootID) then
       t2=MPI_WTime()   

       print *, '  Time to calculate H*X/(H-E)*X:', t2-t1
       print *, '*****************************************'
       print *
    end if
 end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine  progHXFile(x0fname, svMode, x1fname, hxfname)
    character(len=*), intent(IN) :: x0fname,x1fname, hxfname
    logical, intent(IN) :: svMode
    double precision, dimension(plen(sF)) :: X0, X
    double complex,   dimension(plen(sF)) :: X0CX, XCX    
    double precision :: t1, t2
    logical :: log0

    if (id==rootID) then
       print *
       print *, '*****************************************'
       print *, '*        Testing  HX method             *'
       print *, '*****************************************'
       print *       

       print *, ' Read initial X from file:', x0fname

       if (sCX) then
          print *, ' H*X: H is complex, X is complex'
       else
          if (sAP) then
             print *, ' H*X: H is real, X is complex'
          else
             print *, ' H*X: H is real, X is real'
          end if
       end if
       print *, ' Save initial X to file:', x1fname
       print *, ' Save final H*X to file:', hxfname

       select case (sHC)
       case (:TA0)
           print *, '    Calculating H*X'
       case (TA1)
           print *, '    Calculating (H-E)*X'
       case (TAAP)
           print *, '    Calculating (H-E+iAP)*X'
       case (TAAPP)
           print *, '    Calculating (H-E+iAPP)*X'
       case (TAAPR)
           print *, '    Calculating (H-E+iAPR)*X'
       case (TMAP)
           print *, '    Calculating (H-E-iAP)*X'
       case (TMAPP)
           print *, '    Calculating (H-E-iAPP)*X'
       case (TMAPR)
           print *, '    Calculating (H-E-iAPR)*X'
       end select

       t1=MPI_WTime()

    end if

    if (sCX)  then   
         log0 = loadDiag_CX(.true.,x0fname,X0CX)
         log0 = saveDiag_CX(svMode,x1fname,X0CX)
         call HX_CX(plen(sF),X0CX, XCX)
         log0 = saveDiag_CX(svMode,hxfname,XCX)
    else
         if (sAP) then
             log0 = loadDiag_CX(.true.,x0fname,X0CX)
             log0 = saveDiag_CX(svMode,x1fname,X0CX)
             call HX_DX(plen(sF),X0CX, XCX)
             log0 = saveDiag_CX(svMode,hxfname,XCX)
         else
  	     log0 = loadDiag(.true.,x0fname,X0)
             log0 = saveDiag(svMode,x1fname,X0)
             call HX(plen(sF),X0, X)
             log0 = saveDiag(svMode,hxfname,X)
         end if
    end if

    if (id==rootID) then
       t2 = MPI_WTime()  

       print *, '  Time to calculate H*X/(H-E)*X:', t2-t1
       print *, '*****************************************'
    end if
 end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


