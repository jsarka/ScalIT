!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    Get the DVR points for Legendre Polynomial   c
!c                 functions                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
! Input parameters:                                c
!    N : number of DVR nodes                       c
!    X1, X2 : X range                              c
! Output parameters                                c
!    Xout : DVR nodes                              c
!    Wout : DVR weight                             c
!ccccccccccccccccccccccccccccccccccccccccccccccccccc

!*******************************************************
subroutine YjNode(N, Xout)       ! in [-1, 1], cos(theta)
   implicit none
   integer, intent(IN)  :: N          
   double precision, intent(OUT) :: Xout(N)   

   double precision :: beta(N)
   double precision :: vec, work
   integer :: i, ldz, info

   ldz = 1

   DO I = 1, N
       beta(I) =  1.0D0 / DSQRT((2.0D0*i+1.0D0)*(2.0D0*i-1.0D0)) * i 
   END DO

   CALL DSTEV('N', N, XOUT, beta, vec, ldz, work, INFO)

end 

!***********************************************************************
subroutine YjmNode(j, m, Xout)       ! in [-1, 1], cos(theta)
   implicit none
   integer, intent(IN)  :: j, m
   double precision, intent(OUT) :: Xout(J)

   double precision :: beta(J)
   double precision :: vec, work
   integer :: i, ldz, info

   ldz = 1

   DO I = 1, j
       beta(I) =  DSQRT(DBLE((i*i-m*m))/(2.0D0*i+1.0D0)*(2.0D0*i-1.0D0)) 
   END DO

   CALL DSTEV('N', j, XOUT, beta, vec, ldz, work, INFO)

end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine YjNodes(N, x, w) 
   implicit none
   integer, intent(IN) :: N
   double precision, intent(OUT) :: x(N), w(N)

   call Gau_Legendre(-1.0D0, 1.0D0, N, x, w) 

end

subroutine Gau_legendre(x1, x2, n, x, w)
  implicit none
  double precision, parameter  :: EPSI = 1.0D-14
  double precision, parameter  :: PIE  = 3.14159265358979323846264338328d0
  integer, intent(in)          :: n       !number of abscissas and weights
  double precision, intent(in) :: x1, x2  !range

  double precision, intent(out) :: x(N)  !real abscissas 
  double precision, intent(out) :: w(N)  !real weights

  integer :: i,j,m
  double precision :: p1,p2,p3,pp
  double precision :: xl,xm,z,z1
      
  m=(n+1)/2;     xm=0.5d0*(x2+x1);    xl=0.5d0*(x2-x1)

  do i=1,m
     z=cos( PIE * (i-.25d0)/(n+.5d0))
     
     LOOP: do 
        p1=1.d0;         p2=0.d0
        
        do j=1,n
           p3=p2;        p2=p1
           p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
        end do

        pp=n*(z*p1-p2)/(z*z-1.d0)
        z1=z;          z=z1-p1/pp
        
        if(abs(z-z1) <= EPSI )    exit

     end do LOOP
     
     x(i)=xm-xl*z;        x(n+1-i)=xm+xl*z
     w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
     w(n+1-i)=w(i)

  end do

end

!****************************************************************

