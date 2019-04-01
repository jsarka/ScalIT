!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Chebyshev Polynomial, T(n)=2x*T(n-1)-T(n-2), |x|<=1  c 
!c       T0=1,T1=x,T2=2x^2-1.                            c
!c                                                       c
!c Input parameters:                                     c
!c    N: order of Chebeshev Polynomial, N > 1            c
!c    x: the x value, make sure |x|<=1 before call it    c
!c                                                       c
!c Output parameter:                                     c
!c   chev[N+1]:  chev[i+1]=Ti(x)                         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine ChevPolys(N0, N, x, chev)
   implicit none
   integer, intent(IN) :: N0, N
   double precision, intent(IN)  :: x(N0)
   double precision, intent(OUT) :: chev(N0,N+1)

   integer :: i
   double precision :: x2(N0)
   x2(1:N0) = x(1:N0) + x(1:N0)   ! 2x
   chev(1:N0, 1) = 1.0D0
   chev(1:N0, 2) = x(1:N0)

   do i = 3, N+1
      chev(1:N0,i) = x2(1:N0)*chev(1:N0,(i-1))-chev(1:N0,(i-2))
   end do

end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine ChevPoly(N, x, chev)
   implicit none
   integer, intent(IN) :: N
   double precision, intent(IN)  :: x(1) 
   double precision, intent(OUT) :: chev(N+1)

   call ChevPolys(1,N,x,chev)
end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Get the DVR points and weight                  c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine ChevNode(N, x)
   implicit none
   integer, intent(IN) :: N
   double precision, intent(OUT) :: x(N)

   double precision, parameter :: PI = 3.1415926535897932D0
   integer :: i
  
   do i = 1, N
      x(i) = (i-0.5D0) 
   end do
   x(1:N) = dcos(x(1:N)*PI/N)
end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Get the DVR points and weight                  c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine ChevNodes(N, x, w)
   implicit none
   integer, intent(IN) :: N
   double precision, intent(OUT) :: x(N), w(N)

   double precision, parameter :: PI = 3.1415926535897932D0
   integer :: i
  
   do i = 1, N
      x(i) = (i-0.5D0)
   end do
   x(1:N) = dcos(x(1:N)*PI/N)
   w(1:N) = PI / N
end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

