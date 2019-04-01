!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         Subroutine to implement Lanczos algorithm      c
!c                      using Polynomials                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! It returns integer:
!  0: Input parameters error
! -n: not convergent
!  n: successful
!

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function PolyLAN(M0, LAN_EIG)  
     integer, intent(IN)           :: M0    
     double precision, intent(OUT) :: LAN_EIG(M0)
     
     double precision :: X (myLen)
     integer :: LAN   

!cccccccccccccccccccccccc   
     call random_number(X)
     polyLan = LAN(myLen, X, M0, PolyEHX, LAN_EIG)

end function

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function PolyLAN2(M1, LAN_EIG1, M2, LAN_EIG2)                          
     integer, intent(IN)           :: M1, M2                              
     double precision, intent(OUT) :: Lan_Eig1(M1), LAN_EIG2(M2)

     double precision :: X(myLen)
     integer :: LANEIG 
      
     call random_number(X)
     PolyLAN2 = LANEIG(myLen, X, M1, M2, PolyEHX, LAN_EIG1, LAN_EIG2)  

end function 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function PolyLan_CONV( M0, LAN_EIG, RES)
     integer, intent(IN)           :: M0        
     double precision, intent(OUT)  :: LAN_EIG(M0)
     double precision, intent(OUT) :: RES

     integer :: lan_Conv
     integer :: nType, lanStart, lanStep,lanEigStep,lanNMax
     double precision :: X(myLen), lanE0, lanEtol

     nType  = 0; lanE0=sConv%mE0; lanETol=sConv%mTol
     if (sConv%mStart < M0)   then
         lanStart = 2*M0
     else
         lanStart = sConv%mStart
     end if
     if (sConv%mStep  < 1) then
        lanStep = 5
     else
        lanStep = sConv%mStep
     end if
     if (sConv%mGap < 1)   then
        lanEigStep = 1
     else
        lanEigStep = sConv%mGap
     end if
     if (sConv%mMax < lanStart )   then
         lanNMax = lanStart + M0
     else
         lanNMax = sConv%mMax
     end if

     call random_number(X)

     PolyLan_CONV = lan_CONV(lanE0, lanETOL, nType, myLen, X, lanStart, lanStep, &
                   lanEigStep, M0, lanNMax, polyEHX, LAN_EIG, RES)
end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function polyLan_CONVERG(M0, LAN_EIG, RES)            
     integer, intent(IN)           :: M0
     double precision, intent(OUT) :: LAN_EIG(M0)
     double precision, intent(OUT) :: RES
     
     double precision :: X(myLen), lanE0, lanEtol
     integer :: lan_Converg
     integer :: nType, lanStart, lanStep,lanEigStep,lanNMax
     
     nType  = 0; lanE0=sConv%mE0; lanETol=sConv%mTol
     if (sConv%mStart < M0)   then
         lanStart = 2*M0
     else
         lanStart = sConv%mStart
     end if
     if (sConv%mStep  < 1) then
        lanStep = 5
     else
        lanStep = sConv%mStep
     end if
     if (sConv%mGap < 1)   then
        lanEigStep = 1
     else
        lanEigStep = sConv%mGap
     end if
     if (sConv%mMax < lanStart )   then
         lanNMax = lanStart + M0
     else
         lanNMax = sConv%mMax
     end if

     call random_number(X)
                         
     polyLAN_CONVERG = Lan_CONVERG(lanE0, lanETOL, nType, myLen, X, lanStart, &
			  lanStep, lanEigStep, M0, lanNMax, HX, LAN_EIG, RES)  
end function 
!cccccccaccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine setPolyPara(pNum, pCoeff)
    integer, intent(IN) :: pNum 
    double precision, intent(IN) :: pCoeff(pNum)

    if (pNum > MaxPolyNum) then 
          polyNum = MaxPolyNum
          print *, ' PolyNum is greater than the maximum!'
    else
          polyNum = pNum
    end if

    polyCoeff(1:pNum) = pCoeff(1:pNum)

end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PolyEHX(NIN, X, Y)

   integer, intent(IN) :: NIN
   double precision, intent(IN) :: X(NIN)
   double precision, intent(OUT):: Y(NIN)

   integer :: i
   double precision :: tmp1(NIN), tmp2(NIN)

   call EHX(NIN, X, tmp1)   !  tmp1= (H-E0)*X
   Y(1:NIN) = PolyCoeff(1)*tmp1(1:Nin)

   do i = 2, polyNum
      call EHX( NIN, tmp1, tmp2)  ! tmp2=(H-E0)*tmp1
      Y(1:NIN) = Y(1:NIN) + PolyCoeff(i)*tmp2(1:NIN)
      tmp1(1:NIN) = tmp2(1:NIN)
   end do

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
