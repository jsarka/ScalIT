!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Subroutine to the abscissas and weights for       c
!c   Gauss-Chebychev function: See "Numerical Recipes" 	 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function ChevX0(N, X0)  
   implicit none
   integer, intent(IN)  :: N          
   double precision, intent(OUT) :: X0(N)   

   double precision :: beta(N), vec, work
   integer :: i, ldz, info

   ldz = 1;  X0(1:N)=0.0D0
   beta(1)=sqrt(2.0D0)/2.0D0
   do I = 2, N
       beta(I) =  0.5D0
   end do

   call DSTEV('N', N, X0, beta, vec, ldz, work, info)

   chevX0 = (info==0)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function ChevX0W0(N, X0, W0)  
   implicit none
   integer, intent(IN)  :: N          
   double precision, intent(OUT) :: X0(N),W0(N)   
   double precision, parameter  :: PIE  = 3.14159265358979324D0

   double precision :: beta(N), vec(N,N), work(3*N)
   integer :: i, ldz, info

   ldz = N;  X0(1:N)=0.0D0
   beta(1)=sqrt(2.0D0)/2.0D0
   do I = 2, N
       beta(I) =  0.5D0
   end do

   call DSTEV('V', N, X0, beta, vec, ldz, work, info)

   chevX0W0 = (info==0)

   if (ChevX0W0) then
       W0(1:N) = PIE*vec(1,1:N)*vec(1,1:N)
   end if

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function ChevX0Vec(N, X0, Vmat)  
   implicit none
   integer, intent(IN)  :: N          
   double precision, intent(OUT) :: X0(N), Vmat(N, N)

   double precision :: beta(N), work(3*N)
   integer :: i, ldz, info

   ldz = N;  X0(1:N)=0.0D0
   beta(1)=sqrt(2.0D0)/2.0D0
   do I = 2, N
       beta(I) =  0.5D0
   end do

   call DSTEV('V', N, X0, beta, Vmat, ldz, work, info)

   chevX0Vec = (info==0)

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine GauChevX0W0(n, x0, w0)
   implicit none
   integer, intent(in)          :: n   
   double precision, intent(out) :: x0(N), w0(N)

   call Gau_Chev(N, X0, W0)
end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!subroutine Gau_Chev(n, x0, w0)
!   implicit none
!   integer, intent(in)          :: n   
!   double precision, intent(out) :: x0(N), w0(N)
!
!   double precision, parameter  :: PIE  = 3.14159265358979324D0   
!
!   integer :: i
!
!   do I = 1, N
!      x0(I) = DCOS(PIE*(I-0.5D0)/N)
!   end do
!   
!   W0(1:N) = PIE/N
!      
!end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
