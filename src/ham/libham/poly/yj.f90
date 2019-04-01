!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             Normalized Legendre Polynomial,             c
!c   YL=sqrt(2L+1)*(2l-1))*P(L-1)/L                        c
!c             -(1-1/L)sqrt((2l+1)/(2l-3))P(L-2)           c    
!c    |x|<=1   P0=1,P1=x,P2=(3x^2-1)/2                     c
!c                                                         c
!c Input parameters:                                       c
!c    M: dimension of x                                    c
!c    L: order of Legendre Polynomial, L > 1               c
!c    x: the x value, make sure |x|<=1 before call it      c
!c                                                         c
!c Output parameter:                                       c
!c   pj[i][j]:  Pj[i][j]=P(j)(x[i])                        c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine YjPolys(N0, L, x, yj)
   implicit none
   integer, intent(IN) :: N0, L
   double precision, intent(IN)  :: x(N0)
   double precision, intent(OUT) :: yj(N0,L+1)

   double precision, parameter :: SQRT_1_2=0.70710678118655D0
                      ! sqrt(1.0D0/2.0D0)=sqrt(0.5D0)
   double precision, parameter :: SQRT_3_2=1.22474487139159D0
                      ! sqrt(3.0D0/2.0D0)=sqrt(1.5D0)
   integer :: i
   double precision :: re0, re1, re2  

   yj(1:N0, 1) = SQRT_1_2
   yj(1:N0, 2) = SQRT_3_2*x(1:N0)
   do i = 3, L+1 
      re0 = 1.0D0/(i-1.0D0)
      re1 = dsqrt( ( -1.0D0 + i + i ) * ( -3.0D0 + i + i ) )
      re1 = re0 * re1
      re2 = dsqrt( ( -1.0D0 + i + i ) / ( -5.0D0 + i + i ) )
      re2 = re2 * ( 1.0 - re0 )
      yj(1:N0,i) = re1*x(1:N0)*yj(1:N0,(i-1))-re2*yj(1:N0,(i-2))
   end do

end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine yjPoly(L, x, yj)
   implicit none
   integer, intent(IN) :: L
   double precision, intent(IN)  :: x(1) 
   double precision, intent(OUT) :: yj(L+1)

   call yjPolys(1,L,x,yj)
end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc




