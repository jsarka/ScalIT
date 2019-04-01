!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     CRP calculation: See J. Chem. Phys., 96(6), 4412, 1995        c
!c     P = 4*[APr^1/2*(H-iAP-E)^-1*APp*(H+iAP-E)^-1*APr^1/2]         c
!c  Hear the factor 4 is not included                                c
!c         Y = 1/4 P * X                                             c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine OSBCRP(Ncrp, NM, engy, crpVal)
     integer, intent(IN) :: Ncrp, NM
     double precision, intent(IN)  :: engy(Ncrp)
     double precision, intent(OUT) :: crpVal(Ncrp)

     integer :: Lan_DX_MPI, i

     double complex :: X(pLen(sF))
     double precision :: Eig(NM)

     if (allocated(SQ_APR))   deallocate(SQ_APR)
     allocate(SQ_APR(pLen(1)))

     SQ_APR(1:pLen(1)) = DSQRT(APR(1:pLen(1)))
     AP(1:pLen(1)) = APP(1:pLen(1))+APR(1:pLen(1))

     if (sOSB==TOSBW)  then
	 if(.not. initMOSBW()) return
     end if

     crpVal(1:Ncrp)=-1.0D0

     do i = 1, Ncrp
        call randMat_CX(pLen(sF), 1, X)

        if (sCX) then
           if (LAN_DX_MPI(id,rootID, pLen(sF), X, NM, CRP_HXCX, EIG)>0)      &
               crpVal(i) = 4.0D0*sum(EIG)
        else
           if (LAN_DX_MPI(id,rootID, pLen(sF), X, NM, CRP_HXDX, EIG)>0)      &
               crpVal(i) = 4.0D0*sum(EIG)
        end if
     end do

     if (sOSB == TOSBW)  call FinalMOSBW()

     deallocate(SQ_APR)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     P = APr^1/2*[(H+iAP-E)^-1]*APp*[(H-iAP-E)^-1]*APr^1/2         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  CRP_HXDX(NIN, X, Y)   
     integer, intent(IN) :: NIN
     double complex, intent(IN)  :: X(NIN)
     double complex, intent(OUT) :: Y(NIN)

     double precision :: ERR     
     integer :: num
     double complex :: tmp(NIN)

     !cccccccccccccccccccccccccccccccccccccccccccccccc
     !c            Y = APr^(1/2) * X.                c
     !cccccccccccccccccccccccccccccccccccccccccccccccc        
     Y(1:NIN) = SQ_APR(1:NIN)*X(1:NIN)

     !cccccccccccccccccccccccccccccccccccccccccccccccc
     !c tmp = 1/(H-E0-iAP)*Y. nPC can be changed to  c
     !c      using various preconditioners.          c
     !cccccccccccccccccccccccccccccccccccccccccccccccc
     sHC = TMAP;      
     if (sAP) then
         sPC = TMAP ;  
         if (.NOT. createWHCX_DX(sCX)) return
     else
         sPC = TA1 
     end if
     num = OSB_QMRDX(NIN, Y, tmp, ERR) 

     !cccccccccccccccccccccccccccccccccccccccccccccccc
     !c              Y = APp * X.                    c
     !cccccccccccccccccccccccccccccccccccccccccccccccc 
     Y(1:NIN) = APP(1:NIN)*tmp(1:NIN)

     !cccccccccccccccccccccccccccccccccccccccccccccccc
     !c tmp = 1/(H-E0+iAP)*Y. nPC can be changed to  c
     !c      using various preconditioners.          c
     !cccccccccccccccccccccccccccccccccccccccccccccccc
     sHC = TAAP;       
     if (sAP) then
         sPC = TAAP ;  
         if (.NOT. createWHCX_DX(sCX)) return
     else
         sPC = TA1 
     end if  
     num = num + OSB_QMRDX(NIN, Y, tmp, ERR)

     !cccccccccccccccccccccccccccccccccccccccccccccccc
     !c            Y = APr^(1/2) * X.                c
     !cccccccccccccccccccccccccccccccccccccccccccccccc 
     Y(1:NIN) =  SQ_APR(1:NIN)*tmp(1:NIN) 

     if (id==rootID)  &
         print *, ' QMR Iteration for (H+iAP-E)^-1 and (H-iAP-E)^-1:', num

  end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     P = APr^1/2*[(H+iAP-E)^-1]*APp*[(H-iAP-E)^-1]*APr^1/2         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  CRP_HXCX(NIN, X, Y)   
     integer, intent(IN) :: NIN
     double complex, intent(IN)  :: X(NIN)
     double complex, intent(OUT) :: Y(NIN)

     double precision :: ERR     
     integer :: num
     double complex :: tmp(NIN)

     !cccccccccccccccccccccccccccccccccccccccccccccccc
     !c            Y = APr^(1/2) * X.                c
     !cccccccccccccccccccccccccccccccccccccccccccccccc        
     Y(1:NIN) = SQ_APR(1:NIN)*X(1:NIN)

     !cccccccccccccccccccccccccccccccccccccccccccccccc
     !c tmp = 1/(H-E0-iAP)*Y. nPC can be changed to  c
     !c      using various preconditioners.          c
     !cccccccccccccccccccccccccccccccccccccccccccccccc
     sHC = TMAP;  
     if (sAP) then
         sPC = TMAP ;  
         if (.NOT. createWHCX_DX(sCX)) return
     else
         sPC = TA1 
     end if     
     num = OSB_QMRCX(NIN, Y, tmp, ERR) 

     !cccccccccccccccccccccccccccccccccccccccccccccccc
     !c              Y = APp * X.                    c
     !cccccccccccccccccccccccccccccccccccccccccccccccc 
     Y(1:NIN) = APP(1:NIN)*tmp(1:NIN)

     !cccccccccccccccccccccccccccccccccccccccccccccccc
     !c tmp = 1/(H-E0+iAP)*Y. nPC can be changed to  c
     !c      using various preconditioners.          c
     !cccccccccccccccccccccccccccccccccccccccccccccccc
     sHC = TAAP;       
     if (sAP) then
         sPC = TAAP ;  
         if (.NOT. createWHCX_DX(sCX)) return
     else
         sPC = TA1 
     end if

     num = num + OSB_QMRCX(NIN, Y, tmp, ERR)

     !cccccccccccccccccccccccccccccccccccccccccccccccc
     !c            Y = APr^(1/2) * X.                c
     !cccccccccccccccccccccccccccccccccccccccccccccccc 
     Y(1:NIN) =  SQ_APR(1:NIN)*tmp(1:NIN) 

     if (ID==rootID)   &
         print *, ' QMR Iteration for (H+iAP-E)^-1 and (H-iAP-E)^-1:', num

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

