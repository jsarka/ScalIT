!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                       c
!c   Minimization of a function of n variables using Powell Algorithm    c
!c       Modify Powell subroutien from "Numerical Recipes"               c
!c         It will call linmin function.                                 c
!c                                                                       c
!c Input parameters:                                                     c
!c    N : number of variables                                            c
!c    NP: number of directions                                           c
!c   P(NP)  : starting point                                             c
!c   Xi(NP, NP) : initial direction, normally the unit vector            c
!c   ftol   : tolenrence for convergence                                 c
!c   func   : the function                                               c
!c Output parameters:                                                    c
!c   P(NP)  : Final point                                                c
!c   Xi(NP,NP): Final direction                                          c
!c   fret   : the function value at P                                    c
!c The function return the number of iteration,                          c
!c       < 0 means the convergence fails                                 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function powell( N, NP, P, Xi, ftol, func, fret)
      implicit none
      integer, intent(IN) :: N, NP
      double precision, intent(INOUT) :: P(N)       ! starting point
      double precision, intent(INOUT) :: Xi(N, np)  ! initial direction
      double precision, intent(IN)  :: ftol
      double precision, external    :: func
      double precision, intent(OUT) :: fret

      integer, parameter :: NMAX=20, ITMAX=1000
       
      integer          :: iter, i, ibig, j
      double precision :: del, fp, fptt, t, linmin
      double precision, dimension(N)  :: pt, ptt, xit

!***************************

      fret=func(N, p)

      pt(1:N) = p(1:N)

      powell = - ITMAX

      do iter=1, ITMAX 

         fp=fret
         ibig=0
         del=0.0D0
      
         do  i=1,np
             xit(1:N) = xi(1:N, i)
             fptt=fret
        
             fret = linmin(N, p, xit, func)

             if(abs(fptt-fret) > del)then
                 del  = abs(fptt-fret)
                 ibig = i
             endif
         end do

        if(2.0D0*abs(fp-fret) <= ftol*(abs(fp)+abs(fret))) then
             powell = iter
             return
        end if

        ptt(1:N) = 2.0D0*P(1:N) - pt(1:N)
        xit(1:N) = p(1:N) - pt(1:N)
        pt(1:N)  = p(1:N)

        fptt=func(N,ptt)
  
        if(fptt >= fp)   cycle        

        t = 2.0D0*( fp - 2.0D0*fret + fptt ) * ( fp - fret - del)**2 &
            - del * (fp-fptt)**2
   
        if(t >= 0.0D0 )  cycle        

        fret = linmin(N, p, xit, func)
   
        XI(1:N, ibig) = xi(1:N, np)
        Xi(1:N, np) = xit(1:N) 
      end do
  
  end
