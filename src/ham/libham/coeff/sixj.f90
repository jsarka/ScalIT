!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                               c
!c    Subroutines for angular momentum coupling coefficients     c
!c     Algorithms from  SIAM J. Sci. Comput. Vol25, No.4, P1416  c
!c                                                               c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Calculate the sixj symbols                          c
!c       { j1  j2  j3  }      where ji are integers          c
!c       { j4  j5  j6  }                                     c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function isZero6j(j1,j2,j3,j4,j5,j6)
   implicit none
   integer, intent(IN) :: j1, j2, j3, j4, j5, j6

   logical :: isNotTri

   isZero6j =(     isNotTri(j1, j2, j3) .OR. isNotTri(j3,j4,j5) &
              .OR. isNotTri(j1, j5, j6) .OR. isNotTri(j2,j4,j6))

end

!***************************************************************
logical function isNonZero6j(j1,j2,j3,j4,j5,j6)
   implicit none
   integer, intent(IN) :: j1, j2, j3, j4, j5, j6

   logical :: isTri

   isNonZero6j = ( isTri(j1, j2, j3) .AND. isTri(j3,j4,j5) &
             .AND. isTri(j1, j5, j6) .AND. isTri(j2,j4,j6))
end

!**************************************************************
logical function isTri(j1, j2, j3)
   implicit none
   integer, intent(IN) :: j1, j2, j3

   isTri = ( (j1>=ABS(j2-j3)) .AND. (j1 <= (j2+j3)) .AND. & 
             (j2>=ABS(j1-j3)) .AND. (j2 <= (j1+j3)) .AND. &
             (j3>=ABS(j2-j1)) .AND. (j3 <= (j2+j1)) )  
end 

!*************************************************************
logical function isNotTri(j1, j2, j3)
   implicit none
   integer, intent(IN) :: j1, j2, j3

   isNotTri = ( (j1<ABS(j2-j3)) .OR. (j1 > (j2+j3)) .OR. & 
                (j2<ABS(j1-j3)) .OR. (j2 > (j1+j3)) .OR. &
                (j3<ABS(j2-j1)) .OR. (j3 > (j2+j1)) )
end

!**************************************************************
double precision function sixj(j1, j2, j3, j4, j5, j6)
   implicit none
   integer, intent(IN) :: j1, j2, j3, j4, j5, j6

   double precision, allocatable ::  lnn(:)
   integer :: jMax
   logical :: isZero6j
   double precision :: sixjComp

   if (isZero6j(j1, j2, j3, j4, j5, j6)) then
       sixj = 0.0D0
   else
       jmax = Max(j1+j2+j3, j3+j4+j5, j1+j5+j6, j2+j4+j6)
       jmax = jmax + 2
       allocate(lnn(jmax))
       call lnFN(jmax, lnn)
       sixj = sixjComp(j1, j2, j3, j4, j5, j6, jmax, lnn)
       deallocate(lnn)
   end if
end 

!*******************************************************************
double precision function sixjComp(j1, j2, j3, j4, j5, j6, N, lnNF)
    implicit none
    integer, intent(in) :: j1, j2, j3, j4, j5, j6, N
    double precision, intent(IN) :: lnNF(N)

    INTEGER :: j123, j345, j156, j246
    INTEGER :: j1245, j1346, j2356, jmax
    INTEGER :: k, kmin, kmax
    DOUBLE PRECISION :: tmp, factor
   
    sixjComp=0.0D0          

    j123 = j1 + j2 + j3;        j345 = j3 + j4 + j5
    j156 = j1 + j5 + j6;        j246 = j2 + j4 + j6
    j1245 = j1 + j2 + j4 + j5
    j1346 = j1 + j3 + j4 + j6
    j2356 = j2 + j3 + j5 + j6

    kmin = MAX(j123, j156, j246, j345)
    kmax = MIN(j1245, j1346, j2356)
  !  jmax = MAX(j123+1, j345+1, j156+1, j246+1, kmax)
  !  jmax = jmax + 1

    tmp=(lnNF(j1+j2-j3+1)+lnNF(j2+j3-j1+1)+lnNF(j3+j1-j2+1)-lnNF(j123+2)  &
       + lnNF(j3+j4-j5+1)+lnNF(j4+j5-j3+1)+lnNF(j5+j3-j4+1)-lnNF(j345+2)  &
       + lnNF(j1+j5-j6+1)+lnNF(j5+j6-j1+1)+lnNF(j6+j1-j5+1)-lnNF(j156+2)  &
       + lnNF(j2+j4-j6+1)+lnNF(j4+j6-j2+1)+lnNF(j6+j2-j4+1)-lnNF(j246+2)  &
          ) * 0.5D0 
    
    do k = kmin, kmax
        factor =  lnNF(k-j123+1) + lnNF(k-j156+1)  + lnNF(k-j246+1)  &
                + lnNF(k-j345+1) + lnNF(j1245-k+1) + lnNF(j1346-k+1) &
                + lnNF(j2356-k+1)    

        factor = tmp - factor + lnNF(k+2)    ! (k+1)!

        if (k/2*2 == k) then
           sixJComp = sixJComp + EXP(factor)
        else
           sixJComp = sixJComp - EXP(factor)
        end if

    end do

  end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

