!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Associate Legendre Polynomial                 c    
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

subroutine PjmPolys(N0, jmax, m, x, pjm)
   implicit none
   integer, intent(IN) :: N0, jmax, m
   double precision, intent(IN)  :: x(N0)
   double precision, intent(OUT) :: pjm(N0,jmax-m+1)

   integer :: i, ind
   double precision :: fact, x2(1:N0)
    
   x2(1:N0) = sqrt(1.0D0-x(1:N0)*x(1:N0));
   fact = 1.0D0
   pjm(1:N0,1)=1.0D0
   do i = 1, m
      pjm(1:N0, 1) = -fact*pjm(1:N0,1)*x2(1:N0)
      fact = fact + 2.0D0
   end do

   if (jmax==m) return

   pjm(1:N0, 2) = x(1:N0)*(1.0+m+m)*pjm(1:N0, 1)
   if (jmax == (m+1)) return

   do i = m+2, jmax
      ind = i - m + 1
      pjm(1:N0,ind)=((-1.0D0+i+i)*x(1:N0)*pjm(1:N0,ind-1)-    &
                      (i+m-1.0D0)*pjm(1:N0,(ind-2)))/DBLE(i-m)
   end do

end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PjmPoly(jmax, m, x, pjm)
   implicit none
   integer, intent(IN) :: jmax, m
   double precision, intent(IN)  :: x(1) 
   double precision, intent(OUT) :: pjm(jmax-m+1)

   call PjmPolys(1,jmax,m,x,pjm)
end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!             store all Pjm data
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine AllPjmPolys(N0, jmax, X, pjm)
   implicit none
   integer, intent(IN) :: N0, jmax
   double precision, intent(IN)  :: x(N0)
   double precision, intent(OUT) :: pjm(N0, (jmax+1)*(jmax+2)/2)

   integer :: startInd, m

   do m = 0, jmax
      startInd = m*(jmax+jmax+3-m)/2+1
      call PjmPolys(N0, jmax, m, x, pjm(1,startInd))
   end do

end

subroutine AllPjmPoly(jmax, x, pjm)
   implicit none
   integer, intent(IN) :: jmax
   double precision, intent(IN)  :: x(1)
   double precision, intent(OUT) :: pjm((jmax+1)*(jmax+2)/2)

   call AllPjmPolys(1,jmax,x,pjm)

end  
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

