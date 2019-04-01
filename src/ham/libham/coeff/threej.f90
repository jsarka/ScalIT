!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                               c
!c    Subroutines for angular momentum coupling coefficients     c
!c     Algorithms from  SIAM J. Sci. Comput. Vol25, No.4, P1416  c
!c                                                               c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Wigner 3j symbols (j1  j2  j3 )                               c
!c                   (m1, m2, m3 )                               c
!c   The expression is valid only ji, mi are both integers or    c
!c   both half-integers and J=j1+j2+j3 is a integer              c
!c   It is none zero when j1, j2, j3 form a triangle,            c
!c                        m1+m2+m3=0                             c
!c  Calculate threej symbols                                     c
!c      (j1/2,  j2/3,  j3/2)                                     c
!c      (m1/2,  m2/3,  j3/2)                                     c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Calculate the threej symbols                        c
!c       ( j1  j2  j3  )      where ji, mi are integers      c
!c       ( m1  m2  m3  )                                     c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!***********************************************************************
logical function isZero3j(j1,m1,j2,m2,j3,m3)
   implicit none
   integer, intent(IN) :: j1,m1,j2,m2,j3,m3

   isZero3j= (((m1+m2+m3)/=0) .OR. (j1<abs(m1))  .OR.     &
               (j2<abs(m2)) .OR. (j3<abs(m3))    .OR.     &
               (j1>(j2+j3)) .OR. (j1<abs(j2-j3)) .OR.     &
               (j2>(j3+j1)) .OR. (j2<abs(j3-j1)) .OR.     &
               (j3>(j1+j2)) .OR. (j3<abs(j1-j2)) )

end

!***********************************************************************
logical function isNonZero3j(j1,m1,j2,m2,j3,m3)
   implicit none
   integer, intent(IN) :: j1,m1,j2,m2,j3,m3  
 
   isNonZero3j = (((m1+m2+m3)==0) .AND. (j1>=abs(m1))    .AND. &
                   (j2>=abs(m2))  .AND. (j3<abs(m3))     .AND. &
                   (j1<=(j2+j3))  .AND. (j1>=abs(j2-j3)) .AND. &
                   (j2<=(j3+j1))  .AND. (j2>=abs(j3-j1)) .AND. &
                   (j3<=(j1+j2))  .AND. (j3>=abs(j1-j2)) )
end 

!***********************************************************************
double precision function threej(j1, m1, j2, m2, j3, m3)
   implicit none
   integer, intent(IN) :: j1, m1, j2, m2, j3, m3

   double precision :: threejComp, lnn(j1+j2+j3+2)     
   logical :: isZero3j
  
   threej = 0.0D0

   if (isZero3j(j1,m1,j2,m2,j3,m3))  then
       threej = 0.0D0
   else
       call lnFn((j1+j2+j3+2), lnn)
       threej = threejComp(j1,m1,j2,m2,j3,m3,lnn)
   end if

end 

!*********************************************************************
double precision function threejComp(j1,m1,j2,m2,j3,m3, lnNF)
   implicit none
   integer, intent(IN) :: j1, m1, j2, m2, j3, m3
   double precision, intent(IN) :: lnNF(j1+j2+j3+2) 

   integer :: k, kmin, kmax
   integer :: j123, jm1, jm2, j32m1, j31m2, j12m3
   
   double precision :: tmp, factor

   j123  = j1+j2-j3;    jm1   = j1-m1
   jm2   = j2+m2   ;    j32m1 = j3-j2+m1
   j31m2 = j3-j1-m2;    j12m3 = j1-j2-m3

   threeJComp = 0.0D0;
   kmin = max(-j32m1, -j31m2, 0)
   kmax = min(j123, jm1, jm2)

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !c   calculate delta(j1,j2,j3)*sqrt((j1+m1)!(j1-m1)!    c
   !c      *(j2+m2)!*(j2-m2)!(j+m)!*(j-m)!)                c
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   tmp = 0.5D0*(   lnNF(j1+j2-j3+1) + lnNF(j2+j3-j1+1)   & 
                 + lnNF(j3+j1-j2+1) - lnNF(j1+j2+j3+2)   &
                 + lnNF(j1-m1+1)    + lnNF(j1+m1+1)      &
                 + lnNF(j2-m2+1)    + lnNF(j2+m2+1)      &
                 + lnNF(j3-m3+1)    + lnNF(j3+m3+1) )

   do k = kmin, kmax
      factor = lnNF(k+1)     + lnNF(j123-k+1)  + lnNF(jm1-k+1)   &
             + lnNF(jm2-k+1) + lnNF(j32m1+k+1) + lnNF(j31m2+k+1) 

      factor = tmp - factor

      if (k/2*2 == k) then
         threeJComp = threeJComp + exp(factor)
      else
         threeJComp = threeJComp - exp(factor)
      end if

   end do

   if (j12m3/2*2 /= j12m3)    threeJComp = -threeJComp

end
!*******************************************************************
!*****************************************************************





