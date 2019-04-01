!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Create H Matrix in LANCZOS               c
!c           H(I, J) = <VI|H|VJ>                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function lanH0_MPI(MYID, ROOTID, N, X, VJ, M_LOW,   &
                    M_HIGH, HX, ALPHA, BETA)
   implicit none  
   include 'mpif.h'
   integer, intent(IN)   ::  MYID, ROOTID
   integer, intent(IN)   :: N, M_LOW, M_HIGH
   double precision, intent(IN)    :: X(N)
   external             :: HX
   double precision, intent(INOUT) :: VJ(N, M_HIGH+1)
   double precision, intent(INOUT) :: ALPHA(M_HIGH), BETA(M_HIGH)


   double precision :: WJ(N)
   double precision :: TMP
   integer :: J, IERR
   logical :: LAN_GS_ORTH_MPI
   double precision :: NORM_MPI, DOTPROD_MPI
  
   lanH0_MPI= 0

   if ((M_LOW < 1) .or. (M_LOW >= M_HIGH) )  return   

   if (M_LOW == 1) then
      tmp  = NORM_MPI(MPI_COMM_WORLD, N, X(1:N), IERR)  
      if (tmp == 0.0D0)     return      
      VJ(1:N, 1)   = X(1:N)/tmp  
   end if

   LANH0_MPI = M_HIGH

   do J = M_LOW, M_HIGH              ! WJ = H*VJ(1:N, J) 
     call HX(N, VJ(1:N, J), WJ, IERR)

     if (J /= 1)       &
         WJ(1:N)  = WJ(1:N) - BETA(J-1) * VJ(1:N, J-1)
     
     alpha(J) = DOTPROD_MPI(MPI_COMM_WORLD, N, WJ(1:N), VJ(1:N, J), IERR)      
 
     WJ(1:N) = WJ(1:N) - alpha(J) * VJ(1:N, J)

     if (.not. LAN_GS_ORTH_MPI(N, WJ, J, VJ, BETA(J), IERR)) then
         LanH0_MPI = -J;     return
     end if

  end do 

!  if (myid==rootID) print *, 'After loop, LanH0_MPI', lanH0_MPI

end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Create H Matrix in LANCZOS, H(I,J)=<VI|H|VJ>     c
!c    VI is the LANCZOS vectors:                      c          
!c        WJ = (H-E)^-1*VI(I)                         c
!c        VI(I+1) = WJ - Sum((VI(K):WJ)*VI(K))        c
!c                 k = 1, 2, ..., I                   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function lanHij_MPI(MYID, ROOTID, N, X, M, HX, ALPHA, BETA)
   implicit none
   include 'mpif.h'
   integer, intent(IN)   ::  MYID, ROOTID
   integer, intent(IN)   :: N, M
   double precision, intent(IN) :: X(N)
   external             :: HX
   double precision, intent(OUT) :: ALPHA(M), BETA(M)


   double precision ::  VJ(N,M+1), WJ(N), TMP
   integer :: J, IERR
   logical :: LAN_GS_ORTH_MPI  
   double precision :: NORM_MPI, DOTPROD_MPI

   tmp  = NORM_MPI(MPI_COMM_WORLD, N, X, IERR)  !DSQRT(DOT_PRODUCT(X(1:N), X(1:N)))

   lanHIJ_MPI = 0

   if (tmp == 0.0D0)    return

   VJ(1:N, 1)   = X(1:N)/tmp  

   LANHIJ_MPI = M

   do J = 1, M              ! WJ = H*VJ(1:N, J) 

     call HX(N, VJ(1:N, J), WJ, IERR)
 
     if (J /= 1) then
         WJ(1:N) = WJ(1:N) - beta(J-1) * VJ(1:N, J-1)  
     end if
 
     alpha(J) = DOTPROD_MPI(MPI_COMM_WORLD, N, WJ(1:N), VJ(1:N,J), IERR)

     if (.not. LAN_GS_ORTH_MPI(N, WJ, J, VJ, BETA(J), IERR)) then
         LanHij_MPI = -J
         return
     end if     

  end do 

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


