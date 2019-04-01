!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Legendre Polynomial, PL=x*(2-1/L)*P(L-1)-(1-1/L)P(L-2) c 
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

subroutine PjPolys(N0, L, x, pj)
   implicit none
   integer, intent(IN) :: N0, L
   double precision, intent(IN)  :: x(N0)
   double precision, intent(OUT) :: pj(N0,L+1)

   integer :: i
   double precision :: re0, re1, re2
   pj(1:N0, 1) = 1.0D0
   pj(1:N0, 2) = x(1:N0)
  
   do i = 3, L+1 
      re0  = 1.0D0/(i-1); re1 = 2.0D0 - re0; re2 = 1.0D0 - re0
      pj(1:N0,i) = re1*x(1:N0)*pj(1:N0,(i-1))-re2*pj(1:N0,(i-2))
   end do

end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PjPoly(L, x, pj)
   implicit none
   integer, intent(IN) :: L
   double precision, intent(IN)  :: x(1) 
   double precision, intent(OUT) :: pj(L+1)

   call PjPolys(1,L,x,pj)
end






