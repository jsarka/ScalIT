!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                        c
!c         Subroutine to implement PIST algorithm         c
!c                                                        c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! It returns integer:
!  0: Input parameters error
! -n: not convergent
!  n: successful
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function OSB_PIST(M0, PIST_EIG)     
     integer, intent(IN)      :: M0
     double precision, intent(OUT)  :: PIST_EIG(M0)

     integer :: PIST     

     double precision :: X (myLen)

!cccccccccccccccccccccccc        
     call random_number(X)

     OSB_PIST = PIST(myLen, X, M0, HX, OSB_QMR, PIST_EIG)

end function 


!******************************************************************
integer function OSB_PISTCX(M0, PIST_EIG)      
     integer, intent(IN)          :: M0            
     double complex, intent(OUT)  :: PIST_EIG(M0)

     integer :: PIST_CX   
     double complex :: X(myLen) 
!cccccccccccccccccccccccc 

     call randMat_CX(myLen, 1, x)
       
     OSB_PISTCX = PIST_CX(myLen, X, M0, HX_CX, OSB_QMRCX, PIST_EIG)

end function OSB_PISTCX

!************************************************************
