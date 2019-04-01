!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     CRP calculation: See J. Chem. Phys., 96(6), 4412, 1995        c
!c     P = 4*[APr^1/2*(H-iAP-E)^-1*APp*(H+iAP-E)^-1*APr^1/2]         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine osbCRP( Ncrp, NM, engy, crpVal)
    integer, intent(IN)  :: Ncrp, NM
    double precision, intent(IN)  :: engy(Ncrp)
    double precision, intent(OUT) :: crpVal(Ncrp)

    integer :: i, j
    double complex   :: initV(myLen)

    if (sOSB==TOSBW) then
        if (.NOT. initOSBW())  then
             print *, ' Error in OSBW initialization '; return
        end if
    end if

    if (allocated(SQ_APR))  deallocate(SQ_APR)
    allocate(SQ_APR(mylen), stat=i)
    if (i/=0) then 
        print *, '  Error in allocating memory for SQ_APR'
        if (sOSB==TOSBW) call finalOSBW()
        return
    end if

    SQ_APR(1:myLen) = DSQRT(APR(1:myLen))
    AP(1:myLen) = APP(1:myLen)+APR(1:myLen)

    DO I = 1, Ncrp
       call randMat_CX(mylen, 1, initV)   

       sOSBW%mE0 = engy(I)

       crpval(i) = calOSBCRP(mylen, initV, NM) 

       write(*,10) i, engy(i), crpval(i)
       print *
    END DO

    deallocate(SQ_APR)

    if (sOSB==TOSBW) call finalOSBW()

  10  FORMAT('Step:',I5, ' | Energy(eV):', F10.6, ' | CRP:', F15.9)

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  double precision function calOSBCRP(NIN, X, NM)
     integer, intent(IN) :: NIN, NM
     double complex, intent(IN) :: X(NIN)

     double precision :: eig(NM)
     integer :: Lan_DX

     calOSBCRP = 0.0D0     
 
     if (sCX) then
        if (LAN_DX(NIN, X, NM, CRP_HXCX, EIG)>0)      &
          calOSBCRP = 4.0D0*sum(EIG)
     else
        if (LAN_DX(NIN, X, NM, CRP_HXDX, EIG)>0)      &
          calOSBCRP = 4.0D0*sum(EIG)
     end if

  end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


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
         sPC = TMAP;  
         if (.NOT. createWH_DXCX(sCX)) return   
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
         sPC = TAAP; 
         if (.NOT. createWH_DXCX(sCX)) return   
     else
         sPC = TA1 
     end if  
     num = num + OSB_QMRDX(NIN, Y, tmp, ERR)

     !cccccccccccccccccccccccccccccccccccccccccccccccc
     !c            Y = APr^(1/2) * X.                c
     !cccccccccccccccccccccccccccccccccccccccccccccccc 
     Y(1:NIN) =  SQ_APR(1:NIN)*tmp(1:NIN) 

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
         sPC = TMAP; 
         if (.NOT. createWH_DXCX(sCX)) return
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
         sPC = TAAP; 
         if (.NOT. createWH_DXCX(sCX)) return   
     else
         sPC = TA1 
     end if

     num = num + OSB_QMRCX(NIN, Y, tmp, ERR)

     !cccccccccccccccccccccccccccccccccccccccccccccccc
     !c            Y = APr^(1/2) * X.                c
     !cccccccccccccccccccccccccccccccccccccccccccccccc 
     Y(1:NIN) =  SQ_APR(1:NIN)*tmp(1:NIN) 

     print *, ' QMR Iteration for (H+iAP-E)^-1 and (H-iAP-E)^-1:', num

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

