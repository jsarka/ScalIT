!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         Subroutine to implement PIST algorithm         c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   It returns integer:                                  c
!c     0: Input parameters error                          c
!c    -n: not convergent                                  c
!c     n: successful                                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Calculate the bound states                   c
!c           H0 is a real symmetric matrix                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function OSB_PISTCONV(M0, PIST_EIG, RES)   
     integer, intent(IN)          :: M0
     double precision, intent(OUT) :: PIST_EIG (M0)
     double precision, intent(OUT) :: RES

     integer :: PIST_CONV, PISTF_CONV           
     integer :: pistStart, pistStep, pistEigStep,pistNMax
     double precision    :: X(myLen), pistE0, pistETol
!cccccccccccccccccccccccc   

     sHC = TA1; sPC=TA1  
     if (sOSB==TOSBW) then
	 if (.NOT. initOSBW())  then
 	     OSB_PISTCONV = 0; return
         end if
     end if

     call setPISTConvParam(M0, pistE0, pistETol,pistStart, pistStep,    &
                            pistEigStep,pistNMax)

     call random_number(X)
     
     if (sPT==0) then
        OSB_PISTCONV = PIST_CONV(pistE0, pistETOL, sPISTCONVTYPE,    &
			  myLen,X,pistSTART,PISTSTEP,PISTEIGSTEP,    &
                          M0, pistNMax, HX0,OSB_QMR, PIST_EIG,RES)
     else
        OSB_PISTCONV = PISTF_CONV(pistE0, pistETOL, sPISTCONVTYPE,   &
                          myLen,X,pistSTART,PISTSTEP,PISTEIGSTEP,M0, &
                          pistNMax,HX0,OSB_QMR,PIST_EIG,RES,PTFILE)
     end if

     if (sOSB==TOSBW) call finalOSBW()

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Calculate the bound states                   c
!c           H0 is a Hermitian matrix                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function OSB_PISTCONVCX0(M0, PIST_EIG, RES) 
     integer, intent(IN)          :: M0  
     double complex, intent(OUT)  :: PIST_EIG(M0)
     double precision, intent(OUT) :: RES

     integer :: PIST_CONV_CX, PISTF_CONV_CX  
     integer :: pistStart, pistStep, pistEigStep,pistNMax
     double complex   :: X(myLen)
     double precision :: pistE0, pistETol
!cccccccccccccccccccccccc   

     sHC = TA1; sPC=TA1  
     if (sOSB==TOSBW) then
	if (.NOT. initOSBW())  then
 	     OSB_PISTCONVCX0 = 0; return
        end if
     end if

     call setPISTConvParam(M0, pistE0, pistETol,pistStart, pistStep,   &
                            pistEigStep,pistNMax)

     call randMat_CX(mylen, 1, X)

     if (sPT/=0) then
         OSB_PISTCONVCX0=PISTF_CONV_CX(pistE0,pistETOL,sPISTCONVTYPE, &
                     myLen, X, pistSTART, PISTSTEP, PISTEIGSTEP, M0, &
                     pistNMAX, HX0_CX,OSB_QMRCX, PIST_EIG,RES, PTFile)
     else
        OSB_PISTCONVCX0=PIST_CONV_CX(pistE0,pistETOL,sPISTCONVTYPE,  &
                    myLen, X, pistSTART, PISTSTEP, PISTEIGSTEP, M0, &
                    pistNMAX,HX0_CX,OSB_QMRCX, PIST_EIG,RES)
     end if

     if (sOSB==TOSBW) call finalOSBW()

end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Calculate the resonance states:version 1         c
!c           H0 is a real symmetric matrix                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function OSB_PISTCONVDX1(M0, PIST_EIG, RES) 
     integer, intent(IN)          :: M0  
     double complex, intent(OUT)  :: PIST_EIG(M0)
     double precision, intent(OUT) :: RES

     integer :: PIST_CONV_SX, PISTF_CONV_SX , Nt 
     integer :: pistStart, pistStep, pistEigStep,pistNMax
     double complex   :: X(myLen)
     double precision :: pistE0, pistETol
!cccccccccccccccccccccccc   

     sHC = TMAP      ! (H-E0-iAP)*x
     if (sAP) then
        sPC = TMAP   ! (H-E0-iAP)^-1
     else
        sPC = TA1    ! (H-E0)^-1
     end if

 
     if (sOSB==TOSBW) then
	if (.NOT. initOSBW())  then
 	     OSB_PISTCONVDX1 = 0; return
        end if
     end if
     
     call setPISTConvParam(M0, pistE0, pistETol,pistStart, pistStep,   &
                            pistEigStep,pistNMax)

     call randMat_CX(mylen, 1, X)
       
     if (sPT/=0) then
         OSB_PISTCONVDX1=PISTF_CONV_SX(pistE0,pistETOL,sPISTCONVTYPE, &
                     myLen, X, pistSTART, PISTSTEP, PISTEIGSTEP, M0, &
                     pistNMAX, MHX_DX,OSB_QMRDX, PIST_EIG,RES, PTFile)
     else
        OSB_PISTCONVDX1=PIST_CONV_SX(pistE0,pistETOL,sPISTCONVTYPE,  &
                    myLen, X, pistSTART, PISTSTEP, PISTEIGSTEP, M0, &
                    pistNMAX,MHX_DX,OSB_QMRDX, PIST_EIG,RES)
     end if

     if (sOSB==TOSBW) call finalOSBW()

end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      Calculate the resonance states:Version 1          c
!c           H0 is a complex Hermitian matrix             c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function OSB_PISTCONVCX1(M0, PIST_EIG, RES) 
     integer, intent(IN)          :: M0  
     double complex, intent(OUT)  :: PIST_EIG(M0)
     double precision, intent(OUT) :: RES

     integer :: PIST_CONV_SX, PISTF_CONV_SX  
     integer :: pistStart, pistStep, pistEigStep,pistNMax
     double complex   :: X(myLen)
     double precision :: pistE0, pistETol
!cccccccccccccccccccccccc   

     sHC = TMAP      ! (H-E0-iAP)*x
     if (sAP) then
        sPC = TMAP   ! (H-E0-iAP)^-1
     else
        sPC = TA1    ! (H-E0)^-1
     end if

     if (sOSB==TOSBW) then
	if (.NOT. initOSBW())  then
 	     OSB_PISTCONVCX1 = 0; return
        end if
     end if

     call setPISTConvParam(M0, pistE0, pistETol,pistStart, pistStep,   &
                            pistEigStep,pistNMax)

     call randMat_CX(mylen, 1, X)

     if (sPT/=0) then
         OSB_PISTCONVCX1=PISTF_CONV_SX(pistE0,pistETOL,sPISTCONVTYPE, &
                     myLen, X, pistSTART, PISTSTEP, PISTEIGSTEP, M0, &
                     pistNMAX, MHX_CX,OSB_QMRCX, PIST_EIG,RES, PTFile)
     else
        OSB_PISTCONVCX1=PIST_CONV_SX(pistE0,pistETOL,sPISTCONVTYPE,  &
                    myLen, X, pistSTART, PISTSTEP, PISTEIGSTEP, M0, &
                    pistNMAX,MHX_CX,OSB_QMRCX, PIST_EIG,RES)
     end if

     if (sOSB==TOSBW) call finalOSBW()

end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Calculate the resonance states:version 2         c
!c           H0 is a real symmetric matrix                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function OSB_PISTCONVDX2(M0, PIST_EIG, RES) 
     integer, intent(IN)          :: M0  
     double complex, intent(OUT)  :: PIST_EIG(M0)
     double precision, intent(OUT) :: RES

     integer :: PIST_CONV_CX, PISTF_CONV_CX
     integer :: pistStart, pistStep, pistEigStep,pistNMax
     double complex   :: X(myLen)
     double precision :: pistE0, pistETol
!cccccccccccccccccccccccc   

     sHC = TMAP      ! (H-E0-iAP)*x
     if (sAP) then
        sPC = TMAP   ! (H-E0-iAP)^-1
     else
        sPC = TA1    ! (H-E0)^-1
     end if

 
     if (sOSB==TOSBW) then
	if (.NOT. initOSBW())  then
 	     OSB_PISTCONVDX2 = 0; return
        end if
     end if
     
     call setPISTConvParam(M0, pistE0, pistETol,pistStart, pistStep,   &
                            pistEigStep,pistNMax)

     call randMat_CX(mylen, 1, X)
       
     if (sPT/=0) then
         OSB_PISTCONVDX2=PISTF_CONV_CX(pistE0,pistETOL,sPISTCONVTYPE, &
                     myLen, X, pistSTART, PISTSTEP, PISTEIGSTEP, M0, &
                     pistNMAX, MHX_DX,OSB_QMRDX, PIST_EIG,RES, PTFile)
     else
        OSB_PISTCONVDX2=PIST_CONV_CX(pistE0,pistETOL,sPISTCONVTYPE,  &
                    myLen, X, pistSTART, PISTSTEP, PISTEIGSTEP, M0, &
                    pistNMAX,MHX_DX,OSB_QMRDX, PIST_EIG,RES)
     end if

     if (sOSB==TOSBW) call finalOSBW()

end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    Calculate the resonance states: version 2           c
!c           H0 is a complex Hermitian matrix             c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function OSB_PISTCONVCX2(M0, PIST_EIG, RES) 
     integer, intent(IN)          :: M0  
     double complex, intent(OUT)  :: PIST_EIG(M0)
     double precision, intent(OUT) :: RES

     integer :: PIST_CONV_CX, PISTF_CONV_CX  
     integer :: pistStart, pistStep, pistEigStep,pistNMax
     double complex   :: X(myLen)
     double precision :: pistE0, pistETol
!cccccccccccccccccccccccc   

     sHC = TMAP      ! (H-E0-iAP)*x
     if (sAP) then
        sPC = TMAP   ! (H-E0-iAP)^-1
     else
        sPC = TA1    ! (H-E0)^-1
     end if

     if (sOSB==TOSBW) then
	if (.NOT. initOSBW())  then
 	     OSB_PISTCONVCX2 = 0; return
        end if
     end if

     call setPISTConvParam(M0, pistE0, pistETol,pistStart, pistStep,   &
                            pistEigStep,pistNMax)

     call randMat_CX(mylen, 1, X)

     if (sPT/=0) then
         OSB_PISTCONVCX2=PISTF_CONV_CX(pistE0,pistETOL,sPISTCONVTYPE, &
                     myLen, X, pistSTART, PISTSTEP, PISTEIGSTEP, M0, &
                     pistNMAX, MHX_CX,OSB_QMRCX, PIST_EIG,RES, PTFile)
     else
        OSB_PISTCONVCX2=PIST_CONV_CX(pistE0,pistETOL,sPISTCONVTYPE,  &
                    myLen, X, pistSTART, PISTSTEP, PISTEIGSTEP, M0, &
                    pistNMAX,MHX_CX,OSB_QMRCX, PIST_EIG,RES)
     end if

     if (sOSB==TOSBW) call finalOSBW()

end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine setPISTConvParam(M0, pE0, pETol,pStart, pStep, pEigStep,pNMax)
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
