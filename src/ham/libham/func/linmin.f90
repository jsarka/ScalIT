!cccccccccccccccccccccccccccccccccccccccccccccccccc
!c    Subroutine to find the minimal in 1D.       c
!c      See 'Numerical Recipe, Chapt 10'          c
!c Input parameters:                              c
!c   N: Dimensionality                            c
!c   p: initial point                             c
!c   xi: initial direction                        c
!cccccccccccccccccccccccccccccccccccccccccccccccccc

  double precision function linmin(n, p, xi, func)
      implicit none
      integer, intent(in) :: N
      double precision, intent(INOUT) :: p(N), xi(N)
      double precision, external :: func

      double precision, parameter :: TOL = 1.0D-5 
           ! 1.0e-4 for real, 3.0D-8 for double     
      double precision :: ax, bx, xx, xmin
      double precision :: fa, fb, fx, brent

      ax=0.0D0
      xx=1.0D0

      ! Find the 1D range for the minimization 
      call mnbrak(N, p, xi, ax, xx, bx, func, fa, fx, fb)

      ! Find the minimal 
      linmin = brent (N, p, xi, ax, xx, bx, func, TOL, xmin)

      xi(1:N) = xmin * xi(1:N)
      p(1:N)  = p(1:N) + xi(1:N)
 
  end

!*************************************************************************

double precision function brent( N, pi, xi, ax, bx, cx, f, tol, xmin)
      implicit none      
      integer, intent(IN)           :: N
      double precision, intent(IN)  :: pi(N), xi(N)
      double precision, intent(IN)  :: ax, bx, cx, tol
      double precision, external    :: f
      double precision, intent(OUT) :: xmin
    
      integer, parameter :: ITMAX=10000
      double precision, parameter :: CGOLD=0.381966011D0   !(3-sqrt(5))/2
      double precision, parameter :: ZEPS=1.0D-15 

      integer :: iter
      double precision :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
      double precision, dimension(N) :: xtmp

      a=min(ax,cx);       b=max(ax,cx)
      v=bx;      w=v;      x=v;      e=0.0D0

      xtmp(1:N) = pi(1:N) + x * xi(1:N)
      fx=f(N, xtmp)
      fv=fx;      fw=fx

      do iter = 1, ITMAX
        xm=0.5D0*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.0D0*tol1

        if( abs(x-xm) <= (tol2-.5D0*(b-a))) then ! converg
            xmin = x;       brent = fx
            return
        end if   

        if(abs(e) > tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.0D0*(q-r)

          if(q > 0.0D0)    p=-p        

          q=abs(q);     etemp=e;       e=d

          if(abs(p) >= abs(.5D0*q*etemp).or.p <= q*(a-x).or.p >= q*(b-x)) &
               goto 1                   

          d=p/q;             u=x+d
          if((u-a).lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
              goto 2
        endif

1       if(x >= xm) then
          e=a-x
        else
          e=b-x
        endif

        d=CGOLD*e

2       if(abs(d) >= tol1) then
          u=x+d
        else
          u=x+sign(tol1,d)
        endif

        xtmp(1:N) = pi(1:N) + u*xi(1:N)
        fu=f(N, xtmp(1:N))

        if(fu <= fx) then
          if(u >= x) then
            a=x
          else
            b=x
          endif

          v=w;          fv=fw;       w=x
          fw=fx;        x=u;         fx=fu
        else
          if(u < x) then
            a=u
          else
            b=u
          endif

          if(fu <= fw .or. w == x) then
            v=w;         fv=fw
            w=u;         fw=fu
          else if( (fu <= fv) .or. (v == x) .or. (v == w)) then
            v=u;         fv=fu
          endif
        endif
      end do   

      xmin=x;      brent=fx

 end

!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Find the range for later minimization          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mnbrak(N, p, xi, ax, bx, cx, func, fa, fb, fc)
      implicit none
      integer, intent(IN)  :: N
      double precision, intent(IN)    :: p(N), xi(N)
      double precision, intent(INOUT) :: ax, bx, cx
      double precision, intent(OUT)   :: fa, fb, fc
      double precision, external      :: func
      
      double precision, parameter :: GOLD   = 1.618034D0
      double precision, parameter :: GLIMIT = 100.D0
      double precision, parameter :: TINY   = 1.0D-20

      double precision :: dum,fu,q,r,u,ulim
      double precision :: xtmp(N)

      xtmp = p(1:N) + ax * xi(1:N)
      fa=func(N, xtmp)

      xtmp = p(1:N) + bx * xi(1:N)
      fb=func(N, xtmp)

      if(fb > fa)then
        dum=ax;        ax=bx
        bx=dum;        dum=fb
        fb=fa;         fa=dum
      endif

      cx=bx+GOLD*(bx-ax)
      xtmp = p(1:N) + cx * xi(1:N)
      fc=func(N, xtmp)

 1    if (fb >= fc) then
       
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)

        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))

        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx) > 0.0D0) then
          xtmp = p(1:N) + u * xi(1:N)
          fu=func(N, xtmp)

          if(fu < fc)then
            ax=bx;           fa=fb
            bx=u;            fb=fu
            return
          else if(fu > fb)then
            cx=u;            fc=fu
            return
          endif

          u=cx+GOLD*(cx-bx)
          xtmp = p(1:N) + u * xi(1:N)
          fu=func(N, xtmp)

        else if((cx-u)*(u-ulim) > 0.0D0)then

          xtmp = p(1:N) + u * xi(1:N)
          fu=func(N, xtmp)

          if(fu < fc)then

            bx=cx;          cx=u
            u=cx+GOLD*(cx-bx)

            fb=fc;          fc=fu
            xtmp = p(1:N) + u * xi(1:N)
            fu=func(N, xtmp)

          endif

        else if((u-ulim)*(ulim-cx) >= 0.0D0)then
          u=ulim
          xtmp = p(1:N) + u * xi(1:N)
          fu=func(N, xtmp)
        else
          u=cx+GOLD*(cx-bx)
          xtmp = p(1:N) + u * xi(1:N)
          fu=func(N, xtmp)
        endif

        ax=bx;      bx=cx
        cx=u;       fa=fb
        fb=fc;      fc=fu
        goto 1
      end if

  end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

