!
! Minimize of function using simplex, see Numerical Recipe, Chapter 10
!
! Input Parameters:
!   np:   # of parameters
!   ndim: # of parameters to be optimized
!   ftol: convergence criterion
!   P:    intial vertices, it is a (np,ndim+1) matrix
!         each column is one vertex, and
!         we assume that P(ndim+1:np,:) are the same
!
! Output parameters:
!   P:   final vertex
!   iter: # of iterations, <0 when convergence fails
!   
! It returns the minimal values of the final vertex
!

double precision Function simplex(np, ndim, p, ftol, func, iter)
      integer, intent(IN) :: np, ndim
      double precision, intent(IN) :: ftol
      double precision, intent(INOUT) :: P(np, ndim+1)

      double precision, external :: func
      integer,intent(OUT) :: iter

      integer, parameter :: ITMAX=10000  ! 5000
      double precision, parameter :: REF = -1.0D0, ENL=2.0D0, SHK=0.5D0
      double precision, dimension(ndim+1)    :: y

      integer :: ihi, ilo, inhi, mp
      integer :: i, j, m, n
      double precision :: rtol, sum, swap, ysave, ytry
      double precision :: simplexTry
      double precision, dimension(np) :: psum
      
      mp = ndim + 1
      DO I = 1, mp
         Y(I) = func(np, P(1:np, I))
      END DO

      iter=0
1     psum(1:ndim) = p(1:ndim,1)
      do I=2, mp
          psum(1:ndim) = psum(1:ndim) + p(1:ndim,I)
      end do
      psum(ndim+1:np) = p(ndim+1:np, 1)

2     ilo=1
      if (y(1) > y(2)) then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      endif

      do i=1,ndim+1
        if(y(i)<=y(ilo)) ilo=i
        if(y(i) > y(ihi)) then
          inhi=ihi
          ihi=i
        else if(y(i) > y(inhi)) then
          if(i /= ihi) inhi=i
        endif
      end do

      rtol=2.0D0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))

!      print *, 'iter:', iter, ' psum:', psum

      if (rtol <= ftol) then
        swap=y(1)
        y(1)=y(ilo)
        y(ilo)=swap

        do i=1,ndim
          swap=p(i,1)
          p(i,1)=p(i,ilo)
          p(i,ilo)=swap
        end do

        simplex = y(ilo)
        return
      endif

      if (iter >= ITMAX) then
         iter = -iter
         simplex = y(ilo)
         return  
      end if

      iter=iter+2
   
      ytry=simplexTry(np, ndim, p, y, psum, func, ihi, REF)
!      print *,  ' psum1:',psum(1:np)

      if (ytry <= y(ilo)) then
         ytry=simplexTry(np, ndim, p, y, psum, func,ihi,ENL)
!         print *,  ' psum2:',psum(1:np)
      else if (ytry >= y(inhi)) then
        ysave=y(ihi)

        ytry=simplexTry(np, ndim, p, y, psum, func, ihi, SHK)
!        print *,  ' psum3:',psum(1:np)

        if (ytry >= ysave) then
          do  i=1,mp
            if(i /= ilo)then
              psum(1:ndim) = 0.5D0*(p(1:ndim,i)+p(1:ndim,ilo))
              p(1:ndim, i) = psum(1:ndim)
              y(i)=func(np, psum)
            endif
          end do

          iter=iter+ndim
          goto 1

        endif
      else
        iter=iter-1
      endif

      goto 2

      END

!**********************************************************
      double precision FUNCTION simplexTry(np, ndim, p, y, psum, func, ihi, fac)
      INTEGER, intent(IN) :: np, ndim, ihi
      double precision,intent(IN)  :: fac
      double precision, dimension(np, ndim+1),intent(INOUT) :: p
      double precision, dimension(np),intent(INOUT)    :: psum
      double precision, dimension(ndim+1),intent(INOUT)    :: y

      double precision, external :: func

      double precision ::  fac1,fac2,ytry
      double precision, dimension(np) :: ptry

      fac1=(1.0D0-fac)/ndim
      fac2=fac1-fac

      ptry(1:ndim)=psum(1:ndim)*fac1-p(1:ndim,ihi)*fac2
      ptry(ndim+1:np) = psum(ndim+1:np)

      ytry=func(np, ptry)

      if (ytry < y(ihi)) then
        y(ihi)=ytry
        psum(1:ndim)=psum(1:ndim)-p(1:ndim, ihi)+ptry(1:ndim)
        p(1:ndim, ihi)=ptry(1:ndim)
      endif
     
      simplexTry=ytry
     
      END

