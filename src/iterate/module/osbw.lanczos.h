!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                        c
!c         Subroutine to implement Lanczos algorithm      c
!c                                                        c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! It returns integer:
!  0: Input parameters error
! -n: not convergent
!  n: successful
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function LANCZOS(M0, LAN_EIG)  
     integer, intent(IN)           :: M0    
     double precision, intent(OUT) :: LAN_EIG(M0)
     
     double precision :: X(myLen)
     integer :: LAN   

!cccccccccccccccccccccccc   
     call random_number(X)
     LANCZOS = LAN(myLen, X, M0, HX, LAN_EIG)

end function 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function LANCZOS2(M1, LAN_EIG1, M2, LAN_EIG2)                          
     integer, intent(IN)           :: M1, M2                              
     double precision, intent(OUT) :: LAN_EIG1(M1)
     double precision, intent(OUT) :: LAN_EIG2(M2)                 

     double precision  :: X(myLen)
     integer :: LANEIG 
      
     call random_number(X)
     LANCZOS2 = LANEIG(myLen, X, M1, M2, HX, LAN_EIG1, LAN_EIG2)  

end function   

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function LANCZOS_CONV( M0, LAN_EIG, RES)
     integer, intent(IN)           :: M0        
     double precision, intent(OUT) :: LAN_EIG(M0)
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

     LANCZOS_CONV = lan_CONV(lanE0, lanETOL, nType, myLen, X, lanStart, lanStep, &
                   lanEigStep, M0, lanNMax, HX, LAN_EIG, RES)
end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function LANCZOS_CONVERG(M0, LAN_EIG, RES)            
     integer, intent(IN)           :: M0
     double precision, intent(OUT) :: LAN_EIG(M0)
     double precision, intent(OUT) :: RES
     
     double precision :: X(myLen), lanE0, lanETol
     integer :: lan_Converg,lanStart, lanStep,lanEigStep,lanNMax
     integer :: nType
     
     nType    = 0 ; lanE0 = sConv%mE0; lanETol=sConv%mTol
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
                         
     LANCZOS_CONVERG = Lan_CONVERG(lanE0, lanETOL, nType, myLen, X, lanStart, &
			  lanStep, lanEigStep, M0, lanNMax, HX, LAN_EIG, RES)  
end function
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
