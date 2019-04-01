!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                               c
!c    Subroutines for angular momentum coupling coefficients     c
!c     Algorithms from  SIAM J. Sci. Comput. Vol25, No.4, P1416  c
!c                                                               c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Calculate Clesbsch-Gordon Coefficients C(j1m1,j2,m2,j,m)    c
!c   C(j1m1,j2m2|j3m3)=delta(m1+m2,m3)sqrt[(2j3+1)*(s-2j1)!      c
!c       *(s-2j2)!*(2-2j3)!*(j1-m1)!*(j1+m1)!*(j2-m2)!*(j2+m2)!  c
!c       *(j3-m3)!*(j3+m3)!/(s+1)!] *sumk[(-1)^k/{k!*(j1+j2-j3)! c
!c       *(j1-m1-k)!*(j2+m2-k)!*(j3-j2+m1+k)!*(j3-j1-m2+k)!}])   c 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!***********************************************************************
logical function isZeroCG(j1,m1,j2,m2,j3,m3)
   implicit none
   integer, intent(IN) :: j1,m1,j2,m2,j3,m3
   
   logical :: isZero3j
   isZeroCG = isZero3j(j1,m1,j2,m2,j3,-m3)
end

!***********************************************************************
logical function isNonZeroCG(j1,m1,j2,m2,j3,m3)
   implicit none
   integer, intent(IN) :: j1,m1,j2,m2,j3,m3 
 
   logical :: isNonZero3j
   isNonZeroCG = isNonZero3j(j1,m1,j2,m2,j3,-m3) 

end 

!*********************************************************************
double precision function CG(j1, m1, j2, m2, j3, m3)
   implicit none
   integer, intent(IN) :: j1, m1, j2, m2, j3, m3

   double precision    :: threej 
   integer :: j123 

   CG = threej(j1,m1,j2,m2,j3,-m3)*dsqrt(1.0D0+j3+j3)

   j123 = j1-j2+m3
   if (j123/2*2 /= j123) CG = -CG

end 
   
!*********************************************************************
double precision function CGComp(j1, m1, j2, m2, j3, m3, lnNF)
   implicit none
   integer, intent(IN) :: j1, m1, j2, m2, j3, m3

   double precision :: lnNF(j1+j2+j3+2) 

   double precision :: threejComp
   integer :: j123
   
   CGComp = threejComp(j1,m1,j2,m2,j3,-m3,lnNF)*dsqrt(j3+j3+1.0D0)

   j123 = j1-j2+m3
   if (j123/2*2 /= j123) CGComp = -CGComp

end

!**********************************************************************
