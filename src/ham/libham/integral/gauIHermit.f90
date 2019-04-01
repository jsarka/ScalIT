!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Subroutine to the abscissas and weights for       c
!c   Gauss-Hermite function: See "Numerical Recipes" 	 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function HermX0(N, X0)
   implicit none
   integer, intent(IN)  :: N          
   double precision, intent(OUT) :: X0(N)   

   double precision :: beta(N), vec, work
   integer :: i, ldz, info

   ldz = 1;  X0(1:N)=0.0D0
   do I = 1, N
       beta(I) =  DSQRT(0.5D0*I)
   end do

   call DSTEV('N', N, X0, beta, vec, ldz, work, info)
   HermX0 = (info==0)
end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 logical function HermX0W0(N, X0, W0)
   implicit none
   integer, intent(IN)  :: N          
   double precision, intent(OUT):: X0(N),W0(N) 
   double precision, parameter  :: PIM2 = 1.77245385090552D0   ! 1/pi^(1/2)  

   double precision :: beta(N), vec(N,N), work(3*N)
   integer :: i, ldz, info

   ldz = N;  X0(1:N)=0.0D0
   do I = 1, N
       beta(I) =  DSQRT(0.5D0*I)
   end do

   call DSTEV('V', N, X0, beta, vec, ldz, work, info)

   HermX0W0 = (info==0)

   if (HermX0W0)  W0(1:N)=PIM2*vec(1,1:N)*vec(1,1:N)
end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function HermX0Vec(N, X0, Vmat)
   implicit none
   integer, intent(IN)  :: N          
   double precision, intent(OUT):: X0(N), Vmat(N,N) 

   double precision :: beta(N), work(3*N)
   integer :: i, ldz, info

   ldz = N;  X0(1:N)=0.0D0
   do I = 1, N
       beta(I) =  DSQRT(0.5D0*I)
   end do

   call DSTEV('V', N, X0, beta, Vmat, ldz, work, info)

   HermX0Vec = (info==0)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!logical function hermite_alpha(alpha, N, X0)
!   implicit none
!   double precision, intent(IN) :: alpha
!   integer, intent(IN)          :: N          
!   double precision, intent(OUT):: X0(N)
!
!   double precision :: beta(N),vec,work
!   integer :: i, ldz, info
!
!   ldz = 1;  X0(1:N)=0.0D0
!   do I = 1, N
!       beta(I) =  DSQRT(0.5D0*I)
!   end do
!
!   call DSTEV('N', N, X0, beta, vec, ldz, work, INFO)
!
!   hermite_alpha=(info==0)
!   if (hermite_alpha)    X0(1:N) = X0(1:)/alpha  
!end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine GauHermX0W0(N, X0, W0)
   implicit none
   integer, intent(IN)  :: N          
   double precision, intent(OUT):: X0(N),W0(N) 

   call gau_hermite(N, X0, W0)
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!subroutine gau_Hermite(n, x, w)
!   implicit none
!   double precision, parameter  :: EPSI = 1.0D-14
!   double precision, parameter  :: PIM4 = .7511255444649425D0   ! 1/pi^(1/4)
!   integer, parameter :: MAXIT = 100 
!
!   integer, intent(in)          :: n       
!   double precision, intent(out):: x(N), w(N)
!   
!   integer :: i,its,j,m
!   double precision :: p1,p2,p3,pp,z,z1
!
!    m=(n+1)/2
!    do i=1,m
!        if(i == 1)then
!          z = dsqrt(1.0D0*(2*n+1))-1.85575D0/(2*n+1)**(0.16667D0)
!        else if(i == 2)then
!          z=z-1.14*n**.426/z
!        else if (i == 3)then
!          z=1.86*z - .86*x(1)
!        else if (i == 4)then
!          z=1.91*z-.91*x(2)
!        else
!          z=2.*z-x(i-2)
!        endif
!
!        do its=1,MAXIT
!          p1=PIM4;       p2=0.d0
!
!          do j=1,n
!            p3=p2;       p2=p1
!            p1=z*dsqrt(2.d0/j)*p2-dsqrt(dble(j-1)/dble(j))*p3
!          end do
!
!          pp = dsqrt(2.d0*n)*p2
!          z1 = z;            z=z1-p1/pp
!
!          if(abs(z-z1) <= EPSI) exit
!
!        end do
!
!        x(i)=z;               x(n+1-i)=-z
!        w(i)=2.d0/(pp*pp);    w(n+1-i)=w(i)
!
!   end do
!      
!end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
