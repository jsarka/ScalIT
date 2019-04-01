!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Prepare non-orthogonal basis sets for DVR         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         X2=U^T*S0^(-1/2)V^T*X1*V*S0^(-1/2)*U          c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function DVRUX(N, V, S0, U, X1, X2)    
   implicit none
   integer, intent(IN) :: N
   double precision, intent(IN)  :: V(N,N),S0(N),X1(N,N),U(N,N)
   double precision, intent(OUT) :: X2(N,N)

!**************************************************   
   integer :: i, info
   double precision, allocatable  :: X3(:,:)
!**************************************************
   allocate (X3(N,N), stat=info)

   DVRUX = (info==0)

   if (DVRUX) then

      do i = 1, N
         X2(1:N, i) = V(1:N,i)/sqrt(S0(i))
      end do

      call DGEMM('N','N',N,N,N,1.0D0,X2,N,U,N,0.0D0,X3,N)    

      call VTAV(N,N,X1,X3,X2)

      deallocate(X3)

   end if

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         X2=U^H*S0^(-1/2)V^H*X1*V*S0^(-1/2)*U          c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function DVRUXCx(N, VCx, S0,UCx, X1Cx, X2Cx)    
   implicit none
   integer, intent(IN) :: N
   double precision, intent(IN) :: S0(N)
   double complex, intent(IN)   :: VCx(N,N),X1Cx(N,N),UCx(N,N)
   double complex, intent(OUT)  :: X2Cx(N,N)

!**************************************************   
   double complex, parameter :: ONECX=(1.0D0,0.0D0),ZEROCX=(0.0D0,0.0D0)
   integer :: i, info
   double complex, allocatable  :: X3Cx(:,:)
!**************************************************
   allocate (X3Cx(N,N), stat=info)

   DVRUXCx = (info==0)

   if (DVRUXCx) then
      do i = 1, N
         X2Cx(1:N, i) = VCx(1:N,i)/sqrt(S0(i))
      end do

      call ZGEMM('N','N',N,N,N,ONECX,X2Cx,N,UCx,N,ZEROCx,X3Cx,N)      

      call VHAVCx(N,N,X1Cx,X3Cx,X2Cx)
   
      deallocate(X3Cx)

   end if

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         X2=V*S0^(1/2)*U*X1*(V*S0^(1/2)*U)^H           c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function RDVRVX(N, V, S0, U, X1, X2)    
   implicit none   
   integer, intent(IN) :: N
   double precision, intent(IN) :: V(N,N),S0(N),X1(N,N),U(N,N)
   double precision, intent(OUT):: X2(N, N)

!**************************************************   
   integer :: i, info
   double precision, allocatable  :: X3(:,:)
!**************************************************
   allocate (X3(N,N), stat=info)

   RDVRVX = (info==0)

   if (RDVRVX) then
      do i = 1, N
         X2(1:N, i) = V(1:N,i)*sqrt(S0(i))
      end do

      call DGEMM('N','N', N, N, N,1.0D0,X2,N,U,N,0.0D0,X3,N)      

      call VAVT(N, N, X1, X3, X2)
   
      deallocate(X3)
   end if

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         X2=V*S0^(1/2)*U*X1*(V*S0^(1/2)*U)^H           c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function RDVRUXCx(N, VCx, S0, UCx, X1Cx, X2Cx)    
   implicit none   
   integer, intent(IN) :: N
   double precision, intent(IN) :: S0(N)
   double complex, intent(IN)   :: VCx(N,N),X1Cx(N,N),UCx(N,N)
   double complex, intent(OUT)  :: X2Cx(N, N)

!************************************************** 
   double complex, parameter :: ONECX=(1.0D0,0.0D0),ZEROCX=(0.0D0,0.0D0)  
   integer :: i, info
   double complex, allocatable  :: X3Cx(:,:)
!**************************************************
   allocate (X3Cx(N,N), stat=info)

   RDVRUXCx = (info==0)

   if (RDVRUXCx) then
      do i = 1, N
         X2Cx(1:N, i) = VCx(1:N,i)*sqrt(S0(i))
      end do

      call ZGEMM('N','N', N, N, N,ONECX,X2Cx,N,UCx,N,ZEROCx,X3Cx,N)

      call VAVH(N,N,X1Cx,X3Cx,X2Cx)
   
      deallocate(X3Cx)
   end if

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
