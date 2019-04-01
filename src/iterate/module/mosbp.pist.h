!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         Subroutine to implement PIST algorithm         c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   It returns integer:                                  c
!c       0: Input parameters error                        c
!c      -n: not convergent                                c
!c       n: successful                                    c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function OSB_PIST(M0, PIST_EIG)     
     integer, intent(IN)      :: M0
     double precision, intent(OUT)  :: PIST_EIG(M0)

     integer :: PIST_MPI     

     double precision :: X (plen(sF))

!cccccccccccccccccccccccc        
     call random_number(X)

     OSB_PIST = PIST_MPI(pLen(sF),X,M0,HX,OSB_QMR,PIST_EIG)

end function 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function OSB_PISTCX(M0, PIST_EIG)      
     integer, intent(IN)          :: M0            
     double complex, intent(OUT)  :: PIST_EIG(M0)

     integer :: PIST_CX_MPI   
     double complex :: X(plen(sF)) 
!cccccccccccccccccccccccc 

     call randVec_CX(pLen(sF), x)
       
     OSB_PISTCX = PIST_CX_MPI(pLen(sF),X,M0,HX_CX,OSB_QMRCX,PIST_EIG)

end function 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

