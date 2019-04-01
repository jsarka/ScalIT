!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Get the DVR points using Legendre Polynomial  as basis set  c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Input parameters:                                            c
!    N : number of DVR nodes                                   c
!    X1, X2 : X range                                          c
! Output parameters                                            c
!    Xout : DVR nodes                                          c
!    Wout : DVR weight                                         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function LegeX0(N, X0)      
   implicit none
   integer, intent(IN)  :: N          
   double precision, intent(OUT) :: X0(N)   

   double precision :: beta(N),vec,work
   integer :: i, ldz, info

   ldz = 1;  X0(1:N)=0.0D0
   do I = 1, N
       beta(I) =  1.0D0 / DSQRT((2.0D0*i+1.0D0)*(2.0D0*i-1.0D0)) * i 
   end do

   call DSTEV('N', N, X0, beta, vec, ldz, work, info)

   LegeX0 = (info==0) 

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function LegeX0W0(N, X0, W0)       
   implicit none
   integer, intent(IN)  :: N          
   double precision, intent(OUT) :: X0(N), W0(N)   

   double precision :: beta(N), vec(N,N), work(3*N)
   integer :: i, ldz, info

   ldz = N;  X0(1:N)=0.0D0
   do I = 1, N
       beta(I) =  1.0D0 / DSQRT((2.0D0*i+1.0D0)*(2.0D0*i-1.0D0)) * i 
   end do

   call DSTEV('V', N, X0, beta, vec, ldz, work, info)

   LegeX0W0 = (info==0) 

   if (LegeX0W0) then
       W0(1:N) = 2.0D0*vec(1,1:N)*vec(1,1:N)
   end if

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function LegeX0Vec(N, X0, Vmat)       
   implicit none
   integer, intent(IN)  :: N          
   double precision, intent(OUT) :: X0(N), Vmat(N,N)

   double precision :: beta(N), work(3*N)
   integer :: i, ldz, info

   ldz = N;  X0(1:N)=0.0D0
   do I = 1, N
       beta(I) =  1.0D0 / DSQRT((2.0D0*i+1.0D0)*(2.0D0*i-1.0D0)) * i 
   end do

   call DSTEV('V', N, X0, beta, Vmat, ldz, work, info)

   LegeX0Vec = (info==0) 

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine GauLegeX0W0(N, X0, W0) 
   implicit none
   integer, intent(IN) :: N
   double precision, intent(OUT) :: X0(N), W0(N)

   call Gau_Legendre(-1.0D0, 1.0D0, N, X0, W0) 

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!subroutine Gau_legendre(x1, x2, n, x, w)
!  implicit none
!  double precision, parameter  :: EPSI = 1.0D-14
!  double precision, parameter  :: PIE  = 3.14159265358979323846264338328d0
!  integer, intent(in)          :: n       !number of abscissas and weights
!  double precision, intent(in) :: x1, x2  !range
!
!  double precision, intent(out) :: x(N)  !real abscissas 
!  double precision, intent(out) :: w(N)  !real weights
!
!  integer :: i,j,m
!  double precision :: p1,p2,p3,pp
!  double precision :: xl,xm,z,z1
!      
!  m=(n+1)/2;     xm=0.5d0*(x2+x1);    xl=0.5d0*(x2-x1)
!
!  do i=1,m
!     z=cos( PIE * (i-.25d0)/(n+.5d0))
!     
!     LOOP: do 
!        p1=1.d0;         p2=0.d0
!        
!        do j=1,n
!           p3=p2;        p2=p1
!           p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
!        end do
!
!        pp=n*(z*p1-p2)/(z*z-1.d0)
!        z1=z;          z=z1-p1/pp
!        
!        if(abs(z-z1) <= EPSI )    exit
!
!     end do LOOP
!     
!     x(i)=xm-xl*z;        x(n+1-i)=xm+xl*z
!     w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
!     w(n+1-i)=w(i)
!
!  end do
!
!end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!subroutine YjmNode(j, m, Xout)       ! in [-1, 1], cos(theta)
!   implicit none
!   integer, intent(IN)  :: j, m
!   double precision, intent(OUT) :: Xout(J)
!
!   double precision :: beta(J)
!   double precision :: vec, work
!   integer :: i, ldz, info
!
!   ldz = 1
!
!   do I = 1, j
!       beta(I) =  DSQRT(dble((i*i-m*m))/(2.0D0*i+1.0D0)*(2.0D0*i-1.0D0)) 
!   end do
!
!   call DSTEV('N', j, XOUT, beta, vec, ldz, work, INFO)
!
!end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
