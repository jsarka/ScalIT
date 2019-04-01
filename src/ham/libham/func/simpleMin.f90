!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                       c
!c Minimization of a function of n variables without change of dirction  c
!c                                                                       c
!c Input parameters:                                                     c
!c    N : number of variables                                            c
!c    NP: number of directions                                           c
!c   P(N)  : starting point                                              c
!c   Xi(N, NP) : initial direction, normally the unit vector             c
!c   Xr(2, NP) : range at each direction                                 c
!c   ftol   : tolenrence for convergence                                 c
!c   func   : the function                                               c
!c Output parameters:                                                    c
!c   P(N)  : Final point                                                 c
!c   Xi(N,NP): Final direction                                           c
!c   Xr(2,NP): final range                                               c
!c   fret   : the function value at P                                    c
!c The function return the number of iteration,                          c
!c       < 0 means the convergence fails                                 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function simpleMin( N, NP, P, Xi, Xr, ftol, func, fret)
      implicit none
      integer, intent(IN) :: N, NP
      double precision, intent(INOUT) :: P(N)    ! starting point
      double precision, intent(INOUT) :: Xi(N, np)   ! initial direction
      double precision, intent(INOUT)  :: Xr(2,NP)  ! range at that direction
      double precision, intent(IN)  :: ftol
      double precision, external    :: func
      double precision, intent(OUT) :: fret

      integer, parameter ::  ITMAX=1000
       
      integer          :: iter, i 
      double precision :: fpold, fp, lim_linmin
      double precision, dimension(N)  :: pt, ptt, xit

!***************************
      fret=func(N, p)
      simpleMin = 0
      do iter=1, ITMAX 
         fpold = fret
         
         do  i=1, np   ! one iteration for each direction
             xit(1:N) = xi(1:N, i)
        
             fp = lim_linmin(N, p, xit, Xr(1, i), Xr(2,i),func)
                 print *, 'SImple Min, iter=',iter, fp, fp
             if( fp < fret )  fret = fp
             
         end do
         
         simpleMin = simpleMin+1
         print *, 'SImple Min, iter=',iter, fpold, fret, fp

         if (abs(fpold-fret) < ftol)  return
         
      end do

      simpleMin = - simpleMin
  
 end

!cccccccccccccccccccccccccccccccccccccccccccccccccc
!c    Subroutine to find the minimal in 1D.       c
!c      See 'Numerical Recipe, Chapt 10'          c
!c Input parameters:                              c
!c   N: Dimensionality                            c
!c   p: initial point                             c
!c   xi: initial direction                        c
!cccccccccccccccccccccccccccccccccccccccccccccccccc

 double precision function lim_linmin(n, p, xi, x1, x2, func)
      implicit none
      integer, intent(in) :: N
      double precision, intent(INOUT) :: P(N) 
      double precision, intent(OUT)   :: Xi(N)  
      double precision, intent(INOUT) :: x1, x2
      double precision, external   :: func

      double precision, parameter :: TOL = 1.0D-5 
         ! 1.0e-4 for real, 3.0D-8 for double
      integer, parameter :: NMax = 1000           
         ! maximum number to call lim_mnbrak     

      double precision ::  xx, xmin
      double precision ::  brent, lim_mnbrak


      ! Find the 1D range for the minimization 
      xx = lim_mnbrak(N, p, xi, x1, x2, Nmax, func, lim_linmin)

      ! Find the minimal 
      if ((xx == x1) .or. (xx==x2)) then
          xmin = xx          
      else
          lim_linmin = brent (N, p, xi, x1, xx, x2, func, TOL, xmin)
      end if

      p(1:N)  = p(1:N) + xmin*xi(1:N)
 
  end

!*************************************************************************


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Find the minimum range for later brent optimization   c
!c  It is a special case where the range is given, i.e.   c
!c  to find a mid_point, cx, so that ax < cx < bx, and    c
!c  f(c) < f(a), f(b),  random method is used.            c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 double precision function lim_mnbrak(N, p, xi, ax, bx, Nmax, func, fx)
      implicit none
      integer, intent(IN)  :: N, Nmax
      double precision, intent(IN)    :: P(N),Xi(N)
      double precision, intent(INOUT) :: ax, bx
      double precision, external      :: func
      double precision, intent(OUT)   :: fx
      
      double precision :: xtmp(N)
      double precision :: cx, dx, ex, fmin, fa,fb
      integer :: cnt
 
      xtmp(1:N) = p(1:N) + ax * xi(1:N)
      fa=func(N, xtmp)

      xtmp(1:N) = p(1:N) + bx * xi(1:N)
      fb=func(N, xtmp)
 
      fmin=min(fa,fb) 
      cnt = 1;      dx = bx-ax
      call random_seed()
      call random_number(ex)     
      cx = ax + dx*ex

      xtmp(1:N) = p(1:N) + cx * xi(1:N)
      fx=func(N, xtmp)

      do while ((cnt<=Nmax).and.(fx>fmin))
         call random_number(ex)     
         cx = ax + dx*ex

         xtmp(1:N) = p(1:N) + cx * xi(1:N)
         fx=func(N, xtmp)
         cnt = cnt + 1
      end do

      if (fx<fmin) then
         lim_mnbrak = cx
      else
         if (fa<fb) then
            lim_mnbrak = ax
            fx = fa
         else
            lim_mnbrak = bx
            fx = fb
         end if
      end if

      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
