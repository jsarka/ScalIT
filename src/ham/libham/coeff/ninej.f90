!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                               c
!c    Subroutines for angular momentum coupling coefficients     c
!c     Algorithms from  SIAM J. Sci. Comput. Vol25, No.4, P1416  c
!c                                                               c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            calculate a 9-j symbol.                        c
!c   the value of the true arguments in the form             c
!c            { j1 j2 j3 }                                   c
!c            { j4 j5 j6 }   =                               c
!c            { j7 j8 j9 }                                   c
!c                                                           c
!c                       {j1 j4 j7}  {j2 j5 j8} {j3 j6 j9}   c
!c sum(k) [(-1)^2k(2k+1)                                     c
!c                       {j8 j9 k }  {j4 k  j6} {k  j1 j2}   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

logical function isZero9j(a,b,c,d,e,f,g,h,i)
    implicit none
    integer, intent(IN) :: a,b,c,d,e,f,g,h,i 
  
    isZero9j = ( (ABS(a-b)>c) .OR. ((a+b)<c) .OR.  &
                 (ABS(d-e)>f) .OR. ((d+e)<f) .OR.  &
                 (ABS(g-h)>i) .OR. ((g+h)<i) .OR.  &
                 (ABS(a-d)>g) .OR. ((a+d)<g) .OR.  &
                 (ABS(b-e)>h) .OR. ((b+e)<h) .OR.  &
                 (ABS(c-f)>i) .OR. ((c+f)<i) ) 
end  

!*************************************************************
logical function isNonZero9j(a,b,c,d,e,f,g,h,i)
    implicit none
    integer, intent(IN) :: a,b,c,d,e,f,g,h,i
    logical :: isZero9j
 
    isNonZero9j =( .NOT. isZero9j(a,b,c,d,e,f,g,h,i))
end

!***********************************************************
double precision function ninej(a,b,c,d,e,f,g,h,i)
    implicit none
    integer, intent(IN) :: a,b,c,d,e,f,g,h,i

    double precision  :: sixj
    integer :: xlo, xhi,  x
    logical :: isZero9j

    if (isZero9j(a,b,c,d,e,f,g,h,i)) then
        ninej=0.0
    else
        xlo = MAX(ABS(b-f),ABS(a-i),ABS(h-d))
        xhi = MIN(b+f,a+i,h+d)
        ninej = 0.0
        do x = xlo, xhi
           ninej=ninej+(2*x+1)*sixj(a,b,c,f,i,x)*sixj(d,e,f,b,x,h)*  &
                sixj(g,h,i,x,a,d)
        end do
     end if
end

!******************************************************************


