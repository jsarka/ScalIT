!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         Subroutine to implement PIST algorithm         c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   It returns integer:                                  c
!c     0: Input parameters error                          c
!c    -n: not convergent                                  c
!c     n: successful                                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function OSB_PISTCONV(M0, PIST_EIG, dRES)   
     integer, intent(IN)           :: M0
     double precision, intent(OUT) :: PIST_EIG (M0)
     double precision, intent(OUT) :: dRES

     integer :: PIST_CONV_MPI, PISTF_CONV_MPI           
     integer :: pistStart, pistStep, pistEigStep,pistNMax
     double precision    :: X(pLen(sF)), pistE0, pistETol
     double precision    :: PISTConvTime1, PISTConvTime2, DIAGTime1, DIAGTime2
     double precision    :: MPI_Wtime
!cccccccccccccccccccccccc   

     sHC = TA1;     ! H = H0-E0
     sPC = TA1      ! P = (H0-E0)^-1
     OSB_PISTCONV = 0
   
     if (id==rootID)  PISTConvTime1 = MPI_Wtime()
     call setPISTConvParam(M0, pistE0, pistETol,pistStart, pistStep,  &
                            pistEigStep,pistNMax)

     if (sOSB==TOSBW)  then
        if (.NOT. initMOSBW()) return
        if (id==rootID) write(*,100) sOSBW%mCnt, sOSBW%mDe
     end if
     call random_number(X)

     if (id==rootID) then
         PISTConvTime2 = MPI_Wtime()
         print *
         print *, "    ---------------------------------------------------------------------"
         print *, "      MPI Time for PIST: ",PISTConvTime2 - PISTConvTime1
         print *, "    ---------------------------------------------------------------------"
         print *
     end if
     
     if (id==rootID) DIAGTime1 = MPI_Wtime()
     if (sPT==0) then
         OSB_PISTCONV = PIST_CONV_MPI(ID, rootID, pistE0, pistETOL,   &
			   sPISTCONVTYPE, pLen(sF), X, pistSTART,     &
                           PISTSTEP, PISTEIGSTEP, M0, pistNMax, HX0,  &
                           OSB_QMR, PIST_EIG,dRES)
     else
         OSB_PISTCONV = PISTF_CONV_MPI(id, rootid, pistE0, pistETOL,  &
                           sPISTCONVTYPE, pLen(sF), X, pistSTART,     &
                           PISTSTEP,PISTEIGSTEP, M0, pistNMax, HX0,   &
                           OSB_QMR, PIST_EIG, dRES, fPT,              &
                           myRES%gPos(1), myRES%gSize(1))
     end if

    if (sOSB==TOSBW)  call finalMOSBW()
    if (id==rootID) then
         DIAGTime2 = MPI_Wtime()
         print *
         print *, "    ---------------------------------------------------------------------"
         print *, "      MPI Time for Diagonalization: ",DIAGTime2 - DIAGTime1
         print *, "    ---------------------------------------------------------------------"
         print *
    end if
    

    100 format('  The real window size:',I5,' . The adjusted window with:',F10.6)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function OSB_PISTCONVCX0(M0, PIST_EIG, dRES) 
     integer, intent(IN)           :: M0  
     double complex, intent(OUT)   :: PIST_EIG(M0)
     double precision, intent(OUT) :: dRES

     integer :: PIST_CONV_CX_MPI, PISTF_CONV_CX_MPI  
     integer :: pistStart, pistStep, pistEigStep,pistNMax
     double complex   :: X(pLen(sF))
     double precision :: pistE0, pistETol
!cccccccccccccccccccccccc   

     sHC = TA1;     ! H = H0-E0
     sPC = TA1      ! P = (H0-E0)^-1
     OSB_PISTCONVCX0 = 0

     call setPISTConvParam(M0, pistE0, pistETol,pistStart, pistStep,  &
                            pistEigStep,pistNMax)

     if (sOSB==TOSBW) then
        if (.NOT. initMOSBW()) return
        if (id==rootID) write(*,100) sOSBW%mCnt, sOSBW%mDe
     end if

     !call randVec_CX(pLen(sF),X)
     call randMat_CX(plen(SF), 1, X)

     if (sPT/=0) then
        OSB_PISTCONVCX0=PISTF_CONV_CX_MPI(id,rootID,pistE0,pistETOL,  &
             sPISTCONVTYPE,pLen(sF),X,pistSTART,PISTSTEP,PISTEIGSTEP, &
                M0,pistNMAX, HX0_CX,OSB_QMRCX, PIST_EIG,dRES, fPT)
     else
        OSB_PISTCONVCX0=PIST_CONV_CX_MPI(id,rootID,pistE0,pistETOL,   &
             sPISTCONVTYPE,pLen(sF),X,pistSTART,PISTSTEP,PISTEIGSTEP, &
              M0,pistNMAX,HX0_CX,OSB_QMRCX, PIST_EIG,dRES,            &
                           myRES%gPos(1), myRES%gSize(1))
     end if

     if (sOSB==TOSBW)  call finalMOSBW()
     100 format('  The real window size:',I5,' . The adjusted window with:',F10.6)

end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Calculate the Resonance States                     c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function OSB_PISTCONVDX1(M0, PIST_EIG, dRES)   
     integer, intent(IN)           :: M0
     double complex, intent(OUT)   :: PIST_EIG (M0)
     double precision, intent(OUT) :: dRES

     integer :: PIST_CONV_SX_MPI, PISTF_CONV_SX_MPI           
     integer :: pistStart, pistStep, pistEigStep,pistNMax
     double complex   :: X(pLen(sF))
     double precision :: pistE0, pistETol
!cccccccccccccccccccccccc   

     sHC = TMAP;     ! H = H0-E0-iAP
     if (sAP) then
         sPC = TMAP  ! P = (H0-E0-iAP)^-1
     else
         sPC = TA1   ! P = (H0-E0)^-1
     end if
     OSB_PISTCONVDX1 = 0

     call setPISTConvParam(M0, pistE0, pistETol,pistStart, pistStep,    &
                            pistEigStep,pistNMax)

     if (sOSB==TOSBW) then
        if (.NOT. initMOSBW()) return
        if (id==rootID) write(*,100) sOSBW%mCnt, sOSBW%mDe
     end if

     call randMat_CX(plen(SF), 1, X)

     if (sPT==0) then
         OSB_PISTCONVDX1 = PIST_CONV_SX_MPI(ID, rootID, pistE0, pistETOL, &
			   sPISTCONVTYPE, pLen(sF), X, pistSTART,         &
                           PISTSTEP, PISTEIGSTEP, M0, pistNMax, MHX_DX,   &
                           OSB_QMRDX, PIST_EIG,dRES)
     else
         OSB_PISTCONVDX1 = PISTF_CONV_SX_MPI(id, rootid, pistE0, pistETOL,&
                           sPISTCONVTYPE, pLen(sF), X, pistSTART,         &
                           PISTSTEP,PISTEIGSTEP, M0, pistNMax, MHX_DX,    &
                           OSB_QMRDX, PIST_EIG, dRES, fPT,                &
                           myRES%gPos(1), myRES%gSize(1))
     end if

    if (sOSB==TOSBW)  call finalMOSBW()
    100 format('  The real window size:',I5,' . The adjusted window with:',F10.6)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function OSB_PISTCONVDX2(M0, PIST_EIG, dRES)   
     integer, intent(IN)           :: M0
     double complex, intent(OUT)   :: PIST_EIG (M0)
     double precision, intent(OUT) :: dRES

     integer :: PIST_CONV_CX_MPI, PISTF_CONV_CX_MPI           
     integer :: pistStart, pistStep, pistEigStep,pistNMax
     double complex   :: X(pLen(sF))
     double precision :: pistE0, pistETol
     logical :: log0
!cccccccccccccccccccccccc   

     sHC = TMAP;     ! H = H0-E0-iAP
     if (sAP) then
         sPC = TMAP  ! P = (H0-E0-iAP)^-1
     else
         sPC = TA1   ! P = (H0-E0)^-1
     end if
     OSB_PISTCONVDX2 = 0

     call setPISTConvParam(M0, pistE0, pistETol,pistStart, pistStep,    &
                            pistEigStep,pistNMax)

     if (sOSB==TOSBW) then
        if (.NOT. initMOSBW()) return
        if (id==rootID) write(*,100) sOSBW%mCnt, sOSBW%mDe
     end if

     !call random_number(X)
     call randMat_CX(plen(SF), 1, X)

     if (sPT==0) then
         OSB_PISTCONVDX2 = PIST_CONV_CX_MPI(ID, rootID, pistE0, pistETOL, &
			   sPISTCONVTYPE, pLen(sF), X, pistSTART,         &
                           PISTSTEP, PISTEIGSTEP, M0, pistNMax, MHX_DX,   &
                           OSB_QMRDX, PIST_EIG,dRES)
     else
         OSB_PISTCONVDX2 = PISTF_CONV_CX_MPI(id, rootid, pistE0, pistETOL,&
                           sPISTCONVTYPE, pLen(sF), X, pistSTART,         &
                           PISTSTEP,PISTEIGSTEP, M0, pistNMax, MHX_DX,    &
                           OSB_QMRDX, PIST_EIG, dRES, fPT,                &
                           myRES%gPos(1), myRES%gSize(1))
     end if

    if (sOSB==TOSBW)  call finalMOSBW()
    100 format('  The real window size:',I5,' . The adjusted window with:',F10.6)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function OSB_PISTCONVCX1(M0, PIST_EIG, dRES) 
     integer, intent(IN)           :: M0  
     double complex, intent(OUT)   :: PIST_EIG(M0)
     double precision, intent(OUT) :: dRES

     integer :: PIST_CONV_SX_MPI, PISTF_CONV_SX_MPI  
     integer :: pistStart, pistStep, pistEigStep,pistNMax
     double complex   :: X(pLen(sF))
     double precision :: pistE0, pistETol
!cccccccccccccccccccccccc   

     sHC = TMAP;     ! H = H0-E0-iAP
     if (sAP) then
         sPC = TMAP  ! P = (H0-E0-iAP)^-1
     else
         sPC = TA1   ! P = (H0-E0)^-1
     end if
     OSB_PISTCONVCX1 = 0

     call setPISTConvParam(M0, pistE0, pistETol,pistStart, pistStep,  &
                            pistEigStep,pistNMax)

     if (sOSB==TOSBW) then 
       if (.NOT.  initMOSBW()) return
       if (id==rootID) write(*,100) sOSBW%mCnt, sOSBW%mDe
     end if

     call randMat_CX(plen(SF), 1, X)

     if (sPT/=0) then
        OSB_PISTCONVCX1=PISTF_CONV_SX_MPI(id,rootID,pistE0,pistETOL,  &
             sPISTCONVTYPE,pLen(sF),X,pistSTART,PISTSTEP,PISTEIGSTEP, &
                M0,pistNMAX, MHX_CX,OSB_QMRCX, PIST_EIG,dRES, fPT,    &
                           myRES%gPos(1), myRES%gSize(1))
     else
        OSB_PISTCONVCX1=PIST_CONV_SX_MPI(id,rootID,pistE0,pistETOL,   &
             sPISTCONVTYPE,pLen(sF),X,pistSTART,PISTSTEP,PISTEIGSTEP, &
              M0,pistNMAX,MHX_CX,OSB_QMRCX, PIST_EIG,dRES)
     end if

     if (sOSB==TOSBW)  call finalMOSBW()
     100 format('  The real window size:',I5,' . The adjusted window with:',F10.6)

end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function OSB_PISTCONVCX2(M0, PIST_EIG, dRES) 
     integer, intent(IN)           :: M0  
     double complex, intent(OUT)   :: PIST_EIG(M0)
     double precision, intent(OUT) :: dRES

     integer :: PIST_CONV_CX_MPI, PISTF_CONV_CX_MPI  
     integer :: pistStart, pistStep, pistEigStep,pistNMax
     double complex   :: X(pLen(sF))
     double precision :: pistE0, pistETol
!cccccccccccccccccccccccc   

     sHC = TMAP;     ! H = H0-E0-iAP
     if (sAP) then
         sPC = TMAP  ! P = (H0-E0-iAP)^-1
     else
         sPC = TA1   ! P = (H0-E0)^-1
     end if
     OSB_PISTCONVCX2 = 0

     call setPISTConvParam(M0, pistE0, pistETol,pistStart, pistStep,  &
                            pistEigStep,pistNMax)

     if (sOSB==TOSBW) then
        if (.NOT. initMOSBW()) return
        if (id==rootID) write(*,100) sOSBW%mCnt, sOSBW%mDe
     end if

     call randMat_CX(plen(SF), 1, X)
     !call randVec_CX(pLen(sF),X)

     if (sPT/=0) then
        OSB_PISTCONVCX2=PISTF_CONV_CX_MPI(id,rootID,pistE0,pistETOL,  &
             sPISTCONVTYPE,pLen(sF),X,pistSTART,PISTSTEP,PISTEIGSTEP, &
                M0,pistNMAX, MHX_CX,OSB_QMRCX, PIST_EIG,dRES, fPT,    &
                           myRES%gPos(1), myRES%gSize(1))
     else
        OSB_PISTCONVCX2=PIST_CONV_CX_MPI(id,rootID,pistE0,pistETOL,   &
             sPISTCONVTYPE,pLen(sF),X,pistSTART,PISTSTEP,PISTEIGSTEP, &
              M0,pistNMAX,MHX_CX,OSB_QMRCX, PIST_EIG,dRES)
     end if

     if (sOSB==TOSBW)  call finalMOSBW()
     100 format('  The real window size:',I5,' . The adjusted window with:',F10.6)

end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine setPISTConvParam(M0,pE0,pETol,pStart,pStep,pEigStep,pNMax)
     integer, intent(IN) :: M0
     double precision,intent(OUT)   ::  pE0, pETol
     integer,intent(OUT) :: pStart, pStep, pEigStep,pNMax

     pE0 = sConv%mE0; pETol=sConv%mTol

     if (sConv%mStart < M0)   then
         pStart = 2*M0
     else
         pStart = sConv%mStart
     end if

     if (sConv%mStep  < 1) then
        pStep = 5
     else
        pStep = sConv%mStep
     end if

     if (sConv%mGap < 1)   then
        pEigStep = 1
     else
        pEigStep = sConv%mGap
     end if

     if (sConv%mMax < pStart )   then
         pNMax = pStart + M0
     else
         pNMax = sConv%mMax
     end if

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
