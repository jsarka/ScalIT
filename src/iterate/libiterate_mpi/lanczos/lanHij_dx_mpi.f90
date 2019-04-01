!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Create H Matrix in LANCZOS               c
!c           H(I, J) = <VI|H|VJ>                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function lanH0_dx_MPI(MYID, ROOTID, N, X, VJ, M_LOW,   &
                    M_HIGH, HX, ALPHA, BETA)
   implicit none  
   include 'mpif.h'
   integer, intent(IN)   :: MYID, ROOTID
   integer, intent(IN)   :: N, M_LOW, M_HIGH
   double complex, intent(IN)    :: X(N)
   external             :: HX
   double complex,  intent(INOUT) :: VJ(N, M_HIGH+1)
   double precision,intent(INOUT) :: ALPHA(M_HIGH), BETA(M_HIGH)


   double complex   :: WJ(N)
   double precision :: TMP
   integer :: J, IERR
   logical :: LAN_GS_ORTH_cX_MPI
   double precision :: NORM_cx_MPI
   double complex   :: cxtmp, DOTPROD_cx_MPI
  
   lanH0_dx_MPI= 0

   if ((M_LOW < 1) .or. (M_LOW >= M_HIGH) )   return

   if (M_LOW == 1) then
      tmp  = NORM_cx_MPI(MPI_COMM_WORLD, N, X(1:N), IERR)  

      if (tmp == 0.0D0)        return

      VJ(1:N, 1)   = X(1:N)/tmp  
   end if

   if (myid==rootid)  print *, 'before HX:'

   LANH0_dx_MPI = M_HIGH

   do J = M_LOW, M_HIGH              ! WJ = H*VJ(1:N, J) 
     call HX(N, VJ(1:N, J), WJ, IERR)
      
     if (J /= 1)    &
         WJ(1:N)  = WJ(1:N) - BETA(J-1) * VJ(1:N, J-1)

     cxtmp = DOTPROD_cx_MPI(MPI_COMM_WORLD, N, WJ(1:N), VJ(1:N, J), IERR) 
     alpha(J) = dble(cxtmp)  
     WJ(1:N) = WJ(1:N) - alpha(J) * VJ(1:N, J)

     ! Call Full Reorthgonalization
     if (.not. LAN_GS_ORTH_cx_MPI(N, WJ, J, VJ, BETA(J), IERR)) then
         LanH0_dx_MPI = -J
         return
     end if

  end do 

end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Create H Matrix in LANCZOS :  H(I, J) = <VI|H|VJ> c
!c    VI is the LANCZOS vectors:                      c          
!c        WJ = (H-E)^-1*VI(I)                         c
!c        VI(I+1) = WJ - Sum((VI(K):WJ)*VI(K))        c
!c                 k = 1, 2, ..., I                   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function lanHij_dx_MPI( MYID, ROOTID, N, X, M, HX, ALPHA, BETA)
   implicit none
   include 'mpif.h'
   integer, intent(IN)   :: MYID, ROOTID
   integer, intent(IN)   :: N, M
   double complex, intent(IN) :: X(N)
   external             :: HX
   double precision, intent(OUT) :: ALPHA(M), BETA(M)


   double complex ::  VJ(N,M+1),WJ(N)
   double precision :: TMP
   integer :: J, IERR
   logical :: LAN_GS_ORTH_CX_MPI  
   double precision :: NORM_CX_MPI
   double complex   :: cxtmp, DOTPROD_CX_MPI

   tmp  = NORM_CX_MPI(MPI_COMM_WORLD, N, X, IERR)  

   lanHIJ_dx_MPI = 0

   if (tmp == 0.0D0)      return

   VJ(1:N, 1)   = X(1:N)/tmp  

   LANHIJ_dx_MPI = M

   do J = 1, M              ! WJ = H*VJ(1:N, J) 

     call HX(N, VJ(1:N, J), WJ, IERR)
 
     if (J /= 1)            &
         WJ(1:N) = WJ(1:N) - beta(J-1) * VJ(1:N, J-1)     
 
     cxtmp = DOTPROD_CX_MPI(MPI_COMM_WORLD, N, WJ(1:N), VJ(1:N,J), IERR)
     alpha(J) = dble(cxtmp)        

     ! Full Reorthogonization
     if (.not. LAN_GS_ORTH_CX_MPI(N, WJ, J, VJ, BETA(J), IERR)) then
         LanHij_dx_MPI = -J
         return
     end if     

  end do 

end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
