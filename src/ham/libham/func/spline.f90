!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   setup subroutine for cubic spline interpolation   c
!c   Called once before using splint subroutine        c
!c   Adapted from 'Numerical Recipe'                   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine spline(N, x, y, yp1, ypn, y2)
      implicit none
      integer, intent(IN) ::  n
      double precision, intent(IN):: yp1,ypn
      double precision, intent(IN):: X(N),Y(N)
      double precision, intent(OUT) :: Y2(N)
      
      integer :: i,k
      double precision :: p,qn,sig,un
      double precision :: u(N)

      if (yp1 > 0.99D30) then
        y2(1)=0.0D0
        u(1)=0.0D0
      else
        y2(1)=-0.5D0
        u(1)=(3.0D0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.0D0
        y2(i)=(sig-1.0D0)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/    &
             (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do

      if (ypn > 0.99D30) then
        qn=0.0D0
        un=0.0D0
      else
        qn=0.5D0
        un=(3.0D0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0D0)

      do k=n-1,1,-1
          y2(k)=y2(k)*y2(k+1)+u(k)
      end do
     
 end

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
