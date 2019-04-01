!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         Index for CG coefficient                  c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Return # of 1st order indices (j1j2j3m3)      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function getCGJMSize(jmax)
   implicit none
   integer, intent(IN) :: jmax(3)

   integer :: get3jJMSize

   getCGJmSize=get3jJMSize(jmax) 

end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Calculate the total # of CG Coefficients     c
!c       0=<j(i)<=jmax(i)                            c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function getCGCoeffSize(jmax)
   implicit none
   integer, intent(IN) :: jmax(3)

   integer :: get3jCoeffSize
   getCGCoeffSize = get3jCoeffSize(jmax)
end

!************************************************************
integer function getCGCoeffSize1(jmax, jmSize, jLen, jBase)
   implicit none
   integer, intent(IN)  :: jmax(3), jmSize
   integer, intent(OUT) :: jLen(jmSize),jBase(jmSize) 

   integer :: get3jCoeffSize1

   getCGCoeffSize1 = get3jCoeffSize1(jmax,jmSize,jLen,jBase)

end

!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!c          Make sure j1>=j2>=j3>=m3>0              c 
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function getCGPos1(j1,j2,j3,m3)
   implicit none
   integer, intent(IN) :: j1,j2,j3,m3

   getCGPos1 = j1*(6+j1*(11+j1*(6+j1)))/24              &
                 +j2*(2+j2*(3 +j2))/6 + j3*(j3+1)/2 + m3 + 1
end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Calculate # of non-zero CG Coefficients      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getCGLen(jmax, jSize, jLen)
   implicit none
   integer, intent(IN) :: jmax(3), jSize
   integer, intent(OUT):: jLen(jSize)

   call get3jLen(jmax,jSize,jLen)
end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Get (j1j2j3m3) indices for non-zero CG Coefficients     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getCGIndex(jmax, jSize, jInd)
   implicit none
   integer, intent(IN) :: jmax(3), jSize
   integer, intent(OUT):: jInd(jSize)

   call get3jIndex(jmax,jSize,jInd)

end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      Get both the (j1j2j3m3) indices and number of           c
!c               non-zero CG Coefficients                    c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getCGLenInd(jmax, jSize, jLen, jInd)
   implicit none
   integer, intent(IN) :: jmax(3), jSize
   integer, intent(OUT):: jLen(jSize), jInd(4, jSize)

   call get3jLenInd(jmax,jSize,jLen, jInd)

end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c              Get the position of (jm) in the storage             c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function getCGPos( jSize, jLen, jBase, jm)
   implicit none
   integer, intent(IN) :: jSize, jLen(jSize), jBase(jSize), jm(3, 2)

   integer :: rjm(3,2), get3jPos
   
   rjm(1:3,1:2) = jm(1:3,1:2)
   rjm(3,2) = -jm(3,2)
   getCGPos = get3jPos(jSize,jLen,jBase, rjm)

end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Calculate and fill CG Coefficients                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine fillCG(jmax, jmSize, jLen, jBase, jNum, coeff)
   implicit none
   integer, intent(IN) :: jmax(3),jmSize,jNum
   integer, intent(IN) :: jLen(jmSize),jBase(jmSize)
   double precision, intent(OUT) :: coeff(jNum)
 
   double precision :: lnN(jmax(1)+jmax(2)+jmax(3)+2)
   double precision :: cg
   integer :: j1,j2,j3,m1,m2,m3, jm(3,2)
   integer :: j3min, j3max, pos, jsum
   integer :: getCGPos
   
   jsum = sum(jmax(1:3))+2;   call lnFn(jsum, lnN)
   coeff(1:jNum)=0.0D0
   do j1 = 0, jmax(1)
      do j2 = 0, jmax(2)
         j3min=ABS(j1-j2); j3max=min(jmax(3),j1+j2)
         do j3 = j3min, j3max
            do m1 = -j1, j1
               do m2 = -j2, j2
                  m3 = (m1+m2)
                  if (j3 < ABS(m3)) cycle
                    
                  jm(1,1)=j1; jm(2,1)=j2; jm(3,1)=j3
                  jm(1,2)=m1; jm(2,2)=m2; jm(3,2)=m3

                  pos = getCGPos(jmSize, jLen, jBase, jm)

                  if ((pos<1) .OR. (pos > jNum) ) cycle

                  if (coeff(pos) /= 0.0D0)              &
                      coeff(pos)=cg(j1,m1,j2,m2,j3,m3,lnN)
               end do
            end do
         end do
      end do
   end do
end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Reorder the jm(1:3,1:2) where jm[1:3,1]=[j1,j2,j3]     c
!c   and jm[1:3,2]=[m1,m2,m3] so j1>=j2>=j3, and m3>=0      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function ReorderCG(jm, rjm)
   implicit none
   integer, intent(IN)  :: jm(3,2)
   integer, intent(OUT) :: rjm(3,2)

   logical :: Reorder3j
   integer :: rj(3,2)

   rj(1:3,1:2) = jm(1:3,1:2); rj(3,2)=jm(3,2)
   ReorderCG = Reorder3j(rj, rjm)

end 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!*********************************************************
subroutine getCGLen0(jmax, jSize, jLen)
   implicit none
   integer, intent(IN) :: jmax(3), jSize
   integer, intent(OUT):: jLen(jSize)

   call get3jLen0(jmax, jSize, jLen)
end

!*****************************************************************************
subroutine getCGIndex0(jmax, jSize, jInd)
   implicit none
   integer, intent(IN) :: jmax(3), jSize
   integer, intent(OUT):: jInd(4, jSize)

   call get3jIndex0(jmax, jSize, jInd)

end

!*************************************************************************
subroutine getCGLenInd0(jmax, jSize, jLen, jInd)
   implicit none
   integer, intent(IN) :: jmax(3), jSize
   integer, intent(OUT):: jLen(jSize), jInd(4, jSize)

   call get3jLenInd0(jmax, jSize, jLen, jInd)

end

!*******************************************************************************
