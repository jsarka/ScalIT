!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         Subroutine to implement Lanczos algorithm      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function LANCZOS(M0, LAN_EIG)  
     integer, intent(IN)           :: M0    
     double precision, intent(OUT) :: LAN_EIG(M0)
     
     double precision :: X(pLen(sF))
     integer :: LAN_MPI   

!cccccccccccccccccccccccc   
     call random_number(X)
     LANCZOS = LAN_MPI(id, rootID, pLen(sF), X, M0, HX, LAN_EIG)

end function 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function LANCZOS2(M1, LAN_EIG1, M2, LAN_EIG2)                          
     integer, intent(IN)           :: M1, M2                              
     double precision, intent(OUT) :: LAN_EIG1(M1), LAN_EIG2(M2)                 

     double precision  :: X(pLen(sF))
     integer :: LANEIG_MPI 
      
     call random_number(X)
     LANCZOS2 = LANEIG_MPI(id,rootID,pLen(sF),X,M1,M2,HX,LAN_EIG1,LAN_EIG2)  

end function   

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function LANCZOS_CONV( M0, LAN_EIG, dRES)
     integer, intent(IN)           :: M0        
     double precision, intent(OUT) :: LAN_EIG(M0)
     double precision, intent(OUT) :: dRES

     integer :: lan_Conv_MPI
     integer :: nType, lanStart, lanStep,lanEigStep,lanNMax
     double precision :: X(pLen(sF)), lanE0, lanEtol
 
     nType  = 0; 
     call setPISTConvParam(M0, lanE0, lanETol, lanStart, lanStep,     &
                            lanEigStep, lanNMax)

     call random_number(X)

     LANCZOS_CONV=lan_CONV_MPI(id,rootID,lanE0,lanETOL,nType,pLen(sF),&
               X,lanStart,lanStep,lanEigStep,M0,lanNMax,HX,LAN_EIG,dRES)
end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function LANCZOS_CONVERG(M0, LAN_EIG, dRES)            
     integer, intent(IN)           :: M0
     double precision, intent(OUT) :: LAN_EIG(M0)
     double precision, intent(OUT) :: dRES
     
     double precision :: X(pLen(sF)), lanE0, lanETol
     integer :: lan_Converg_MPI,lanStart, lanStep,lanEigStep,lanNMax
     integer :: nType
     
     nType    = 0 
     call setPISTConvParam(M0, lanE0, lanETol, lanStart, lanStep,     &
                            lanEigStep, lanNMax)

     call random_number(X)
                         
     LANCZOS_CONVERG=Lan_CONVERG_MPI(id,rootID,lanE0,lanETOL,nType,   &
          pLen(sF),X,lanStart,lanStep,lanEigStep,M0,lanNMax,HX,LAN_EIG,dRES)
end function
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
