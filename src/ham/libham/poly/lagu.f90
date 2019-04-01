!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Laguerre Polynomial, nL(n)=(2n-1-x)*L(n-1)-(n-1)L(n-2)c 
!c       L0=1,L1=-x+1,L2=(x^2-4x+2)/2.                   c
!c                                                       c
!c Input parameters:                                     c
!c    N: order of Chebeshev Polynomial, N > 1            c
!c    x: the x value, make sure x>=0 before call it      c
!c                                                       c
!c Output parameter:                                     c
!c   lagu[N+1]:  lagu[i+1]=Ln(x)                         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine LaguPolys(N0, N, x, lagu)
   implicit none
   integer, intent(IN) :: N0, N
   double precision, intent(IN)  :: x(N0)
   double precision, intent(OUT) :: lagu(N0,N+1)

   integer :: i
   
   lagu(1:N0, 1) = 1.0D0
   lagu(1:N0, 2) = 1.0D0 - x(1:N0)

   do i = 3, N+1
      lagu(1:N0,i) = (i+i-3-x(1:N0))*lagu(1:N0,(i-1))-(i-2)*lagu(1:N0,(i-2))
      lagu(1:N0,i) = lagu(1:N0,i)/(i-1)
   end do

end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine laguPoly(N, x, lagu)
   implicit none
   integer, intent(IN) :: N
   double precision, intent(IN)  :: x(1) 
   double precision, intent(OUT) :: lagu(N+1)

   call laguPolys(1,N,x,lagu)
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Get DVR points and weight                   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine laguNode(N, x)
   implicit none
   integer, intent(IN)  :: N
   double precision, intent(OUT) :: X(N)

   double precision :: beta(N)
   double precision :: vec, work
   integer :: i, ldz, info

   ldz = 1

   do I = 1, N
       x(I) = 2.0D0*I -1.0D0
       beta(I) = -I
   end do

   call DSTEV('N', N, x, beta, vec, ldz, work, INFO)
end

subroutine laguNodes(N, x, w)
   implicit none
   integer, intent(IN)  :: N
   double precision, intent(OUT) :: X(N), w(N)

   call gau_Laguerre(0.0D0, N, x, w)

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


