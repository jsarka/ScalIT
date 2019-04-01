!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         Subroutine to implement Lanczos algorithm      c
!c                      using Polynomials                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function PolyLAN(M0, LAN_EIG)  
     integer, intent(IN)           :: M0    
     double precision, intent(OUT) :: LAN_EIG(M0)
     
     double precision :: X (pLen(sF))
     integer :: LAN_MPI   

!cccccccccccccccccccccccc   
     call random_number(X)
     polyLan = LAN_MPI(id,rootID,pLen(sF),X,M0,PolyEHX,LAN_EIG)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function PolyLAN2(M1, LAN_EIG1, M2, LAN_EIG2)                          
     integer, intent(IN)           :: M1, M2                              
     double precision, intent(OUT) :: Lan_Eig1(M1), LAN_EIG2(M2)

     double precision :: X(pLen(sF))
     integer :: LANEIG_MPI 
      
     call random_number(X)
     PolyLAN2 = LANEIG_MPI(id,rootID,pLen(sF),X,M1,M2,PolyEHX,  &
                           LAN_EIG1,LAN_EIG2)  

end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function PolyLan_CONV( M0, LAN_EIG, dRES)
     integer, intent(IN)           :: M0        
     double precision, intent(OUT) :: LAN_EIG(M0)
     double precision, intent(OUT) :: dRES

     integer :: lan_Conv_MPI
     integer :: nType, lanStart, lanStep,lanEigStep,lanNMax
     double precision :: X(pLen(sF)), lanE0, lanEtol

     nType  = 0
     call setPISTConvParam(M0, lanE0, lanETol, lanStart, lanStep,      &
                            lanEigStep, lanNMax)

     call random_number(X)

     PolyLan_CONV=lan_CONV_MPI(id,rootid,lanE0,lanETOL,nType, pLen(sF),&
        X,lanStart,lanStep,lanEigStep,M0,lanNMax,polyEHX,LAN_EIG,dRES)

end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function polyLan_CONVERG(M0, LAN_EIG, dRES)            
     integer, intent(IN)           :: M0
     double precision, intent(OUT) :: LAN_EIG(M0)
     double precision, intent(OUT) :: dRES
     
     double precision :: X(pLen(sF)), lanE0, lanEtol
     integer :: lan_Converg_MPI
     integer :: nType, lanStart, lanStep,lanEigStep,lanNMax
     
     nType  = 0
     call setPISTConvParam(M0, lanE0, lanETol, lanStart, lanStep,     &
                            lanEigStep, lanNMax)

     polyLAN_CONVERG=Lan_CONVERG_MPI(id,rootid,lanE0,lanETOL,nType,   &
            pLen(sF),X,lanStart,lanStep,lanEigStep,M0,lanNMax,HX,     &
            LAN_EIG,dRES)  

end function 
!cccccccaccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccaccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
!cccccccaccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PolyEHX(N, X, Y)
   integer, intent(IN) :: N
   double precision, intent(IN) :: X(N)
   double precision, intent(OUT):: Y(N)

   integer :: i
   double precision :: tmp1(N), tmp2(N)

   call EHX(N, X, tmp1)   !  tmp1= (H-E0)*X
   Y(1:N) = PolyCoeff(1)*tmp1(1:N)

   do i = 2, polyNum
      call EHX( N, tmp1, tmp2)  ! tmp2=(H-E0)*tmp1
      Y(1:N) = Y(1:N) + PolyCoeff(i)*tmp2(1:N)
      tmp1(1:N) = tmp2(1:N)
   end do

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
