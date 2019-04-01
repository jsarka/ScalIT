!cccccccccccccccccccccccccccccccccccccccccccccccccc
!c  A subroutine to find the root of a function   c
!cccccccccccccccccccccccccccccccccccccccccccccccccc

double precision function getRoot(X0, X1, func, epsi, success)
    implicit none
    double precision, intent(IN) ::X0, X1, epsi
    double precision, external :: func
    logical, intent(OUT) :: success

    double precision :: x, xstart, xend
    double precision :: fx, fx0, fx1
 
    success = .true.
    xstart = Min(X0, X1)
    xend   = Max(X0, X1)
    fx0 = Func(xstart)
    fx1 = Func(xend)

    DO while ((xend-xstart) > epsi) 
       x  = (xstart+xend)/2.0D0
       fx = func(x)
       if (fx*fx0 < 0) THEN
          xend = x
          fx1  = fx
       else 
          if (fx*fx1<0) THEN
             xstart = x
             fx0 = fx
          else
              success = .false.
              exit
          end if
       end if 
    END DO

    getRoot = (xstart+xend)/2.0D0

end function getRoot
!ccccccccccccccccccccccccccccccccccccccccccccccccc
