!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine splint(N, xa, ya, y2a, M, x, y)
      implicit none
      integer, intent(IN) :: N,M
      double precision, intent(IN) :: xa(n),y2a(n),ya(n)
      double precision, intent(IN) :: X(M)
      double precision, intent(OUT):: Y(M)

      integer :: i, k,khi,klo
      double precision  :: a,b,h

      klo=1;      khi=n

      do I = 1, M
        klo=1;        khi=n
        do while ((khi-klo) > 1) 
           k=(khi+klo)/2
           if(xa(k) > x(I))then
              khi=k
           else
              klo=k
           endif
        end do
      
        h=xa(khi)-xa(klo)
        if (h == 0.0D0 )  return
        a=(xa(khi)-x(I))/h
        b=(x(I)-xa(klo))/h
        y(I)=a*ya(klo)+b*ya(khi)+((a*a*a-a)*y2a(klo)+   &
             (b*b*b-b)*y2a(khi))*(h*h)/6.0D0
      end do
end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
