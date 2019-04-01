!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Create H Matrix in LANCZOS, H(I, J) = <VI|H|VJ>   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function lanH0_DX(N, X, VJ, M_LOW, M_HIGH, HXCX, ALPHA, BETA)
   implicit none
   integer, intent(IN)   :: N, M_LOW, M_HIGH
   double complex, intent(IN)    :: X(N)
   external             :: HXCX
   double complex, intent(INOUT) :: VJ(N, M_HIGH+1)
   double precision,intent(INOUT):: ALPHA(M_HIGH), BETA(M_HIGH)

   double complex  :: WJ(N), CXTMP
   double precision :: TMP
   integer :: J
   logical :: LAN_MGS_ORTH_CX
  
   lanH0_DX = 0

   if ((M_LOW < 1) .or. (M_LOW >= M_HIGH) )     return

   if (M_LOW == 1) then
      tmp  = DSQRT(dble(dot_product(X(1:N), X(1:N))))
      if (tmp == 0.0D0)        return      
      VJ(1:N, 1)   = X(1:N)/tmp  
   end if

   LANH0_DX = M_HIGH

   do J = M_LOW, M_HIGH              ! WJ = H*VJ(1:N, J) 
     call HXCX(N, VJ(1:N, J), WJ)
      
     if (J /= 1)  WJ(1:N)  = WJ(1:N) - BETA(J-1) * VJ(1:N, J-1)     
     
     cxtmp = dot_product(WJ(1:N), VJ(1:N,J))
     alpha(J) = dble(cxtmp)

     WJ(1:N) = WJ(1:N) - alpha(J) * VJ(1:N, J)

     ! Call Full Reorthgonalization
     if (.not. LAN_MGS_ORTH_CX(N, WJ, J, VJ, BETA(J))) then
         LanH0_DX = -J;         return
     end if

  end do 

end function LANH0_DX



!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Create H Matrix in LANCZOS               c
!c       Complex version for Hermitian matrix         c
!c H(I, J) = <VI|H|VJ>, VI is the LANCZOS vectors:    c          
!c        WJ = (H-E)^-1*VI(I)                         c
!c        VI(I+1) = WJ - Sum((VI(K):WJ)*VI(K))        c
!c                 k = 1, 2, ..., I                   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function lanHij_dx(N, X, M, HXcx, ALPHA, BETA)
   implicit none
   integer, intent(IN)   :: N, M
   double complex, intent(IN) :: X(N)
   external             :: HXCX
   double precision,intent(OUT) :: ALPHA(M), BETA(M)

   double complex :: VJ(N, M+1),WJ(N), CXTMP
   double precision :: TMP   
   integer :: J
   logical :: LAN_MGS_ORTH_CX

   tmp  = DSQRT(dble(dot_product(X(1:N), X(1:N))))

   lanHIJ_dx = 0

   if (tmp == 0.0D0)     return

   VJ(1:N, 1)   = X(1:N)/tmp  

   LANHIJ_dx = M

   do J = 1, M              ! WJ = H*VJ(1:N, J) 

     call HXCX(N, VJ(1:N, J), WJ)
 
     if (J /= 1)  WJ(1:N) = WJ(1:N) - beta(J-1) * VJ(1:N, J-1)  
 
     CXTMP = dot_product(WJ(1:N), VJ(1:N,J))
     alpha(J) = dble(CXTMP)

     ! Full Reorthogonization
     if (.not. LAN_MGS_ORTH_CX(N, WJ, J, VJ, BETA(J))) then
         LanHij_DX = -J;          return
     end if     

  end do 

end function LANHIJ_DX



