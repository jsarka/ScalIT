!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Create H Matrix in LANCZOS                  c
!c           H(I, J) = <VI|H|VJ>                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function lanH0(N, X, VJ, M_LOW, M_HIGH, HX, ALPHA, BETA)
   implicit none
   integer, intent(IN)   :: N, M_LOW, M_HIGH
   double precision, intent(IN)   :: X(N)
   external             :: HX
   double precision, intent(INOUT) :: VJ(N, M_HIGH+1)
   double precision, intent(INOUT) :: ALPHA(M_HIGH), BETA(M_HIGH)

   double precision :: WJ(N), TMP
   integer :: J
   logical :: LAN_MGS_ORTH
  
   lanH0= 0

   if ((M_LOW < 1) .or. (M_LOW >= M_HIGH) )   return

   if (M_LOW == 1) then
      tmp  = DSQRT(dot_product(X(1:N), X(1:N)))
      if (tmp == 0.0D0)    return   
      VJ(1:N, 1)   = X(1:N)/tmp  
   end if

   LANH0 = M_HIGH

   do J = M_LOW, M_HIGH              ! WJ = H*VJ(1:N, J) 
     call HX(N, VJ(1:N, J), WJ)
      
     if (J /= 1)   WJ(1:N)  = WJ(1:N) - BETA(J-1) * VJ(1:N, J-1)
     
     alpha(J) = dot_product(WJ(1:N), VJ(1:N,J))
     WJ(1:N) = WJ(1:N) - alpha(J) * VJ(1:N, J)

     ! Call Full Reorthgonalization
     if (.not. LAN_MGS_ORTH(N, WJ, J, VJ, BETA(J))) then
         LanH0 = -J;         return
     end if

  end do 

end function LANH0


!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Create H Matrix in LANCZOS                   c
!c  H(I, J) = <VI|H|VJ> , VI is the LANCZOS vectors:  c          
!c        WJ = (H-E)^-1*VI(I)                         c
!c        VI(I+1) = WJ - Sum((VI(K):WJ)*VI(K))        c
!c                 k = 1, 2, ..., I                   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function lanHij(N, X, M, HX, ALPHA, BETA)
   implicit none
   integer, intent(IN)   :: N, M
   double precision, intent(IN) :: X(N)
   external             :: HX
   double precision, intent(OUT) :: ALPHA(M), BETA(M)

   double precision :: VJ(N, M+1), WJ(N), TMP
   integer :: J
   logical :: LAN_MGS_ORTH  

   tmp  = DSQRT(dot_product(X(1:N), X(1:N)))

   lanHIJ = 0

   if (tmp == 0.0D0)  return

   VJ(1:N, 1)   = X(1:N)/tmp  

   LANHIJ = M

   do J = 1, M              ! WJ = H*VJ(1:N, J) 

     call HX(N, VJ(1:N, J), WJ)
 
     if (J /= 1)  WJ(1:N) = WJ(1:N) - beta(J-1) * VJ(1:N, J-1)     
 
     alpha(J) = dot_product(WJ(1:N), VJ(1:N,J))

     ! Full Reorthogonization
     if (.not. LAN_MGS_ORTH(N, WJ, J, VJ, BETA(J))) then
         LanHij = -J;          return
     end if     

  end do 

end function LANHIJ

