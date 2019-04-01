!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                     c
!c     Subroutine to the abscissas and weights for     c
!c   Gauss-Laguerre function: See "Numerical Recipes"  c
!c                                                     c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Get the DVR points for Laguerre Polynimial        c
!Ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function LaguX0(N, x0)
   implicit none
   integer, intent(IN)  :: N
   double precision, intent(OUT) :: X0(N)

   double precision :: beta(N),vec, work
   integer :: i, ldz, info

   ldz = 1;  X0(1:N)=0.0D0
   DO I = 1, N
       x0(I) = 2.0D0*I -1.0D0
       beta(I) = -DBLE(I)
   END DO

   CALL DSTEV('N', N, X0, beta, vec, ldz, work, info)

   LaguX0 = (info==0)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function LaguX0W0(N, x0, w0)
   implicit none
   integer, intent(IN)  :: N
   double precision, intent(OUT) :: X0(N),W0(N)

   double precision :: beta(N), vec(N,N), work(3*N)
   integer :: i, ldz, info

   ldz = N;  X0(1:N)=0.0D0
   DO I = 1, N
       x0(I) = 2.0D0*I -1.0D0
       beta(I) = -DBLE(I)
   END DO

   CALL DSTEV('V', N, X0, beta, vec, ldz, work, info)

   LaguX0W0 = (info==0)

   if (LaguX0W0) then
       W0(1:N)=vec(1,1:N)*vec(1,1:N)
   end if
end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function LaguX0Vec(N, x0, Vmat)
   implicit none
   integer, intent(IN)  :: N
   double precision, intent(OUT) :: X0(N), Vmat(N,N)

   double precision :: beta(N), work(3*N)
   integer :: i, ldz, info

   ldz = N;  X0(1:N)=0.0D0
   DO I = 1, N
       x0(I) = 2.0D0*I -1.0D0
       beta(I) = -DBLE(I)
   END DO

   CALL DSTEV('V', N, X0, beta, Vmat, ldz, work, info)

   LaguX0Vec = (info==0)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine GauLaguX0W0(N,X0, W0)
   implicit none
   integer, intent(IN)  :: N
   double precision, intent(OUT) :: X0(N),W0(N)

   call gau_Laguerre(0.0D0,N,X0,W0)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!SUBROUTINE gau_Laguerre(alf, n, x, w)
!  implicit none
!  double precision, parameter  :: EPSI = 1.0D-14
!  integer, parameter :: MAXIT = 400 
!  
!  double precision, intent(in) :: alf   
!  integer, intent(in)          :: n     
!  double precision, intent(out):: X(n)  
!  double precision, intent(out):: W(n)  
!  
!  integer :: i,its,j
!  double precision  :: gammln, ai
!  double precision  :: p1,p2,p3,pp,z,z1
!
!!cccccccccccccccccccccccccc
!  do i=1,n
!     
!     if(i == 1)then
!        z=(1.+alf)*(3.+.92*alf)/(1.+2.4*n+1.8*alf)
!     else if(i == 2)then
!        z=z+(15.+6.25*alf)/(1.+.9*alf+2.5*n)
!     else
!        ai=i-2
!        z=z + ((1.+2.55*ai)/(1.9*ai)+1.26*ai*alf/(1.+3.5*ai))   &
!           *(z-x(i-2))/(1.+.3*alf)
!     endif
!
!     do its=1,MAXIT        
!        p1=1.d0;        p2=0.d0
!        
!        do j=1,n
!           p3=p2;       p2=p1
!           p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j
!        end do
!        
!        pp=(n*p1-(n+alf)*p2)/z
!        z1=z;    z=z1-p1/pp
!        
!        if(dabs(z-z1)<= EPSI) exit
!        
!     end do
!     
!     x(i)=z
!     w(i)=-exp(gammln(alf+n)-gammln(1.0D0*n))/(pp*n*p2)
!  end do
!      
!end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccc
!c       Calculate ln(Gama(x))		     c
!ccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!double precision function GAMMLN(xx)   ! ln(gamma(xx))
!  implicit none  
!  double precision, dimension(6), parameter :: cof =   &
!       (/ 76.18009172947146d0, -86.50532032941677d0,   &
!       24.01409824083091d0, -1.231739572450155d0,      &
!       .1208650973866179d-2,-.5395239384953d-5/)  
!  double precision, parameter :: stp =2.5066282746310005d0
!  double precision, intent(in) :: xx  
!
!  integer :: j
!  double precision :: ser,tmp,x,y
!
!!cccccccccccccccccccccccccccc  
!  x=xx;        y=x
!  tmp=x+5.5d0
!  tmp=(x+0.5d0)*log(tmp)-tmp
!  ser=1.000000000190015d0
!    
!  do j=1,6
!     y=y+1.d0
!     ser=ser+cof(j)/y
!  end do
!
!  gammln=tmp+log(stp*ser/x)
!    
!end
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
