!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Normalized Associate Legendre Polynomial          c
!c  (j-m)P(j,m)=x*(2j-1)*P(j-1,m)-(j+m-1)P(j-2,m)          c 
!c    |x|<=1   P(m,m)=(-1)^m(2m-1)!!(1-x2)^(m/2)           c
!c             p(m+1,m)=x(2m+1)P(m,m)                      c
!c                                                         c
!c Input parameters:                                       c
!c    N0: dimension of x                                   c
!c    jmax: order of Legendre Polynomial, jmax>=m>=0       c
!c    x: the x value, make sure |x|<=1 before call it      c
!c                                                         c
!c Output parameter:                                       c
!c   pj[i][j]:  Pj[i][j]=P(j)(x[i])                        c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine YjmPolys(N0, jmax, m, x, yjm)
   implicit none
   integer, intent(IN) :: N0, jmax, m
   double precision, intent(IN)  :: x(N0)
   double precision, intent(OUT) :: yjm(N0,jmax-m+1)

   double precision, parameter :: SQRT_1_2 = 0.70710678118655   ! sqrt(1/2)
   integer :: i, ind
   double precision :: fact, j1, j2, jm1, jm2, x2(1:N0)

   x2(1:N0) = sqrt(1.0D0-x(1:N0)*x(1:N0))
   yjm(1:N0, 1) = SQRT_1_2 
   fact = 1.0D0
   do i = 1, m
      fact = fact + 2.0   ! 2m+1
      yjm(1:N0, 1) = -sqrt(fact/(fact-1.0D0))*yjm(1:N0,1)*x2(1:N0)
   end do

   if (jmax==m) return

   yjm(1:N0, 2) = sqrt(3.0D0+m+m)*x(1:N0)*yjm(1:N0, 1)
   if (jmax == (m+1)) return
   
   do i = m+2, jmax
      ind = i - m + 1
      fact = DBLE(i+i)  ! 2j
      jm1  = DBLE((i-m)*(i+m))
      jm2  = (i-m-1.0D0)*(i+m-1.0D0)
      yjm(1:N0,ind) = sqrt((fact*fact-1.0D0)/jm1)*x(1:N0)*yjm(1:N0,ind-1) &
                 - sqrt(((fact+1.0D0)/(fact-3.0D0))*(jm2/jm1))*yjm(1:N0,(ind-2))
   end do

end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine YjmPoly(jmax, m, x, yjm)
   implicit none
   integer, intent(IN) :: jmax, m
   double precision, intent(IN)  :: x(1) 
   double precision, intent(OUT) :: yjm(jmax-m+1)

   call YjmPolys(1,jmax,m,x,yjm)
end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!             store all Yjm data
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine AllYjmPolys(N0, jmax, X, yjm)
   implicit none
   integer, intent(IN) :: N0, jmax
   double precision, intent(IN)  :: x(N0)
   double precision, intent(OUT) :: yjm(N0, (jmax+1)*(jmax+2)/2)

   integer :: startInd, m

   do m = 0, jmax
      startInd = m*(jmax+jmax+3-m)/2+1
      call YjmPolys(N0, jmax, m, x, yjm(1,startInd))
   end do

end

subroutine AllYjmPoly(jmax, x, yjm)
   implicit none
   integer, intent(IN) :: jmax
   double precision, intent(IN)  :: x(1)
   double precision, intent(OUT) :: yjm((jmax+1)*(jmax+2)/2)

   call AllYjmPolys(1,jmax,x,yjm)

end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


