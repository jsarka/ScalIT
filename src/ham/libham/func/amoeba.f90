!ccccccccccccccccccccccccccccccccccccccccccccccc
!c     Minimize of function using simplex      c     
!c       see Numerical Recipe, Ch 10           c
!ccccccccccccccccccccccccccccccccccccccccccccccc

subroutine amoeba(p,y,mp,np,ndim,ftol,funk,iter)
      integer, intent(IN) :: mp,ndim,np
      double precision, intent(IN) :: ftol
      double precision, intent(INOUT) :: P(mp,np), Y(mp) 
      double precision, external :: funk
      integer,intent(OUT) :: iter

      integer, parameter :: NMAX=20,ITMAX=10000

      integer :: i,ihi,ilo,inhi,j,m,n
      double precision :: rtol, sum, swap, ysave, ytry, amotry
      double precision :: psum(NMAX) 

      iter=0
1     do  n=1,ndim
        sum=0.0D0
        do m=1,ndim+1
          sum=sum+p(m,n)
        end do
        psum(n)=sum
     end do

2     ilo=1
      if (y(1) > y(2)) then
        ihi=1;    inhi=2
      else
        ihi=2;    inhi=1
      endif

      do i=1,ndim+1
        if(y(i)<=y(ilo)) ilo=i
        if(y(i) > y(ihi)) then
          inhi=ihi;     ihi=i
        else if(y(i) > y(inhi)) then
          if(i /= ihi) inhi=i
        endif
      end do
 
      rtol=2.0D0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))

      if (rtol.lt.ftol) then
        swap=y(1);    y(1)=y(ilo)
        y(ilo)=swap

        do n=1,ndim
          swap=p(1,n);  p(1,n)=p(ilo,n)
          p(ilo,n)=swap
        end do

        return
      endif

      if (iter >= ITMAX) return  

      iter=iter+2
   
      ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,-1.0D0)

      if (ytry <= y(ilo)) then
          ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2.0D0)
      else if (ytry >= y(inhi)) then
        ysave=y(ihi)
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5D0)

        if (ytry >= ysave) then
          do  i=1,ndim+1
            if(i /= ilo)then
              do j=1,ndim
                psum(j)=0.5*(p(i,j)+p(ilo,j))
                p(i,j)=psum(j)
              end do
              y(i)=funk(psum)
            endif
          end do

          iter=iter+ndim
          goto 1

        endif
      else
        iter=iter-1
      endif

      goto 2

  end

!**********************************************************
  double precision function amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)
      integer, intent(IN) :: ihi,mp,ndim,np
      double precision,intent(IN)  :: fac
      double precision, intent(INOUT) :: p(mp,np), pSum(mp),y(mp)

      double precision, external :: funk
      integer, parameter :: NMAX = 20

      integer j
      double precision ::  fac1,fac2,ytry
      double precision :: ptry(NMAX)

      fac1=(1.0D0-fac)/ndim
      fac2=fac1-fac

      ptry(1:ndim)=psum(1:ndim)*fac1-p(ihi,1:ndim)*fac2

      ytry=funk(ptry)

      if (ytry < y(ihi)) then
        y(ihi)=ytry
        psum(1:ndim)=psum(1:ndim)-p(ihi,1:ndim)+ptry(1:ndim)
        p(ihi,1:ndim)=ptry(1:ndim)
      endif
     
      amotry=ytry
     
  end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

