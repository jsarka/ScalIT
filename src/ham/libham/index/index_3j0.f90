!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         Index for Gaunt coefficient               c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Return # of 1st order indices (j1j2j3m3)      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function get3jJMSize(jmax)
   implicit none
   integer, intent(IN) :: jmax(3)

   integer :: kmax, kmed, kmin  
   integer :: get3jPos1

   call mmm(jmax(1),jmax(2),jmax(3), kmax, kmed, kmin)

   get3jJmSize=get3jPos1(kmax,kmed,kmin,kmin) 

end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Calculate the total # of 3j Coefficients     c
!c       0=<j(i)<=jmax(i)                            c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function get3jCoeffSize(jmax)
   implicit none
   integer, intent(IN) :: jmax(3)

   integer :: jSize, get3jJMSize
   integer,allocatable :: jLen(:) 

   jSize = get3jJMSize(jmax)

   allocate(jLen(jSize))

   call get3jLen(jmax,jSize,jLen)

   get3jCoeffSize = SUM(jLen)

   deallocate (jLen)
end

!************************************************************
integer function get3jCoeffSize1(jmax, jmSize, jLen, jBase)
   implicit none
   integer, intent(IN)  :: jmax(3), jmSize
   integer, intent(OUT) :: jLen(jmSize),jBase(jmSize) 

   integer :: i

   call get3jLen(jmax,jmSize,jLen)
   jBase(1) = 1
   do i = 1, jmSize-1
      jBase(i+1) = jBase(i)+jLen(i)
   end do

   get3jCoeffSize1 = jBase(jmSize)+jLen(jmSize)-1
end

!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!c          Make sure j1>=j2>=j3>=m3>0              c 
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function get3jPos1(j1,j2,j3,m3)
   implicit none
   integer, intent(IN) :: j1,j2,j3,m3

   get3jPos1 = j1*(6+j1*(11+j1*(6+j1)))/24              &
                 +j2*(2+j2*(3 +j2))/6 + j3*(j3+1)/2 + m3 + 1
end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Calculate # of non-zero 3j Coefficients      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine get3jLen(jmax, jSize, jLen)
   implicit none
   integer, intent(IN) :: jmax(3), jSize
   integer, intent(OUT):: jLen(jSize)

   integer :: j1, j2, j3, m3
   integer :: kmax, pos
   logical :: isTri

   kmax = max(jmax(1),jmax(2),jmax(3))
   pos = 1 ; jLen(1:jSize)=0 
   do j1=0, kmax
      do j2 = 0, j1
         do j3 = 0, j2
            do m3 = 0, j3
               if (isTri(j1,j2,j3))  jLen(pos) = j2+min(j2,j1-m3)+1
               pos = pos + 1
               if (pos > jSize) return
            end do
         end do
      end do
   end do
end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Get (j1j2j3m3) indices for non-zero 3j Coefficients     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine get3jIndex(jmax, jSize, jInd)
   implicit none
   integer, intent(IN) :: jmax(3), jSize
   integer, intent(OUT):: jInd(4, jSize)

   integer :: j1, j2, j3, m3
   integer :: kmax, pos

   kmax = max(jmax(1),jmax(2),jmax(3))
   pos = 1
   do j1=0, kmax
      do j2 = 0, j1
         do j3 = 0, j2
            do m3 = 0, j3
               jInd(1,pos) = j1; jInd(2,pos) = j2
               jInd(3,pos) = j3; jInd(4,pos) = m3
               pos = pos + 1
               if (pos > jSize) return
            end do
         end do
      end do
   end do  
end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      Get both the (j1j2j3m3) indices and number of           c
!c               non-zero 3j Coefficients                    c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine get3jLenInd(jmax, jSize, jLen, jInd)
   implicit none
   integer, intent(IN) :: jmax(3), jSize
   integer, intent(OUT):: jLen(jSize), jInd(4, jSize)

   integer :: j1, j2, j3, m3
   integer :: kmax, pos, jsum
   logical :: isTri

   kmax = max(jmax(1),jmax(2),jmax(3))
   pos = 1; jLen(1:jSize)=0
   do j1=0, kmax
      do j2 = 0, j1
         do j3 = 0, j2
            jsum = j1 + j2 + j3
            do m3 = 0, j3
               jInd(1,pos) = j1; jInd(2,pos) = j2
               jInd(3,pos) = j3; jInd(4,pos) = m3
               if (isTri(j1,j2,j3) ) jLen(pos)=j2+min(j2,j1-m3)+1
               pos = pos + 1
               if (pos > jSize) return
            end do
         end do
      end do
   end do  
end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c              Get the position of (jm) in the storage             c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function get3jPos( jSize, jLen, jBase, jm)
   implicit none
   integer, intent(IN) :: jSize, jLen(jSize), jBase(jSize), jm(3, 2)

   integer :: rjm(3,2), pos, m2m, get3jPos1
   logical :: isZero3j, Reorder3j, ro

   get3jPos = 0;   pos = sum(jm(1:3,1));  m2m = sum(jm(1:3,2))
   if (isZero3j(jm(1,1),jm(1,2),jm(2,1),jm(2,2),jm(3,1),jm(3,2))) &
       return

   ro = Reorder3j(jm, rjm)

   m2m = min(rjm(2,1), rjm(1,1)-rjm(3,2))

   pos = get3jPos1(rjm(1,1), rjm(2,1),rjm(3,1),rjm(3,2))

   if ( (pos>jSize) .OR. (pos<1) .OR. (jLen(pos)==0) .OR.    &
        (rjm(2,2) > m2m) )     return

   get3jPos = jbase(pos) + rjm(2,2) + rjm(2,1) + 1
   if (.NOT. ro) get3jPos = -get3jPos

end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Calculate and fill 3j Coefficients                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine fill3j(jmax, jmSize, jLen, jBase, jNum, coeff)
   implicit none
   integer, intent(IN) :: jmax(3),jmSize,jNum
   integer, intent(IN) :: jLen(jmSize),jBase(jmSize)
   double precision, intent(OUT) :: coeff(jNum)
 
   double precision :: lnN(jmax(1)+jmax(2)+jmax(3)+2)
   double precision :: threej
   integer :: j1,j2,j3,m1,m2,m3, jm(3,2)
   integer :: j3min, j3max, pos, jsum
   integer :: get3jPos
   
   jsum = sum(jmax(1:3))+2;   call lnFn(jsum, lnN)
   coeff(1:jNum)=0.0D0
   do j1 = 0, jmax(1)
      do j2 = 0, jmax(2)
         j3min=ABS(j1-j2); j3max=min(jmax(3),j1+j2)
         do j3 = j3min, j3max
            do m1 = -j1, j1
               do m2 = -j2, j2
                  m3 = -(m1+m2)
                  if (j3 < ABS(m3)) cycle
                    
                  jm(1,1)=j1; jm(2,1)=j2; jm(3,1)=j3
                  jm(1,2)=m1; jm(2,2)=m2; jm(3,2)=m3

                  pos = get3jPos(jmSize, jLen, jBase, jm)

                  if ((pos<1) .OR. (pos > jNum) ) cycle

                  if (coeff(pos) /= 0.0D0)              &
                      coeff(pos)=threej(j1,m1,j2,m2,j3,m3,lnN)
                           
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
logical function Reorder3j(jm, rjm)
   implicit none
   integer, intent(IN)  :: jm(3,2)
   integer, intent(OUT) :: rjm(3,2)

   integer :: jmax, jmin, jmed, jsum

   jsum = sum(jm(1:3,1))
   if (jm(1,1) >= jm(2,1)) then
       jmax = 1; jmin = 2
       Reorder3j = .TRUE.
   else
       jmax = 2; jmin = 1
       Reorder3j = .FALSE.
   end if

   if (jm(jmax,1) > jm(3, 1) ) then
      if (jm(3, 1) < jm(jmin, 1)) then
          jmed = jmin; jmin = 3
      else
          jmed = 3
          Reorder3j = (.NOT. Reorder3j)
      end if
   else
      jmed = jmax; jmax = 3
   end if 
  
   rjm(1,1:2) = jm(jmax,1:2)
   rjm(2,1:2) = jm(jmed,1:2)
   rjm(3,1:2) = jm(jmin,1:2)

   if (rjm(3,2) < 0 )  then
       rjm(1:3,2) = -rjm(1:3,2)
       if (jSum/2*2 /= jSum) Reorder3j = (.NOT. Reorder3j)
   end if
end 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       Reorder (j1,j2,j3) so (k1>=k2>=k3)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!*********************************************************
subroutine get3jLen0(jmax, jSize, jLen)
   implicit none
   integer, intent(IN) :: jmax(3), jSize
   integer, intent(OUT):: jLen(jSize)

   integer :: j1, j2, j3, m1, m2, m3, pos, get3jPos1
   integer :: j3max, j3min, jsum, jm(3,2), rjm(3,2)

   logical :: ltemp, Reorder3j

   jLen(1:jSize)=0
   do j1=0, jmax(1)
      do j2 = 0, jmax(2)
         j3min = ABS(j1-j2); j3max = min(jmax(3), j1+j2)
         do j3 = j3min, j3max
            jsum = j1 + j2 + j3
            if (jsum/2*2 /= jsum) cycle         
            do m1 = -j1, j1 
               do m2 = -j2, j2
                  m3 = -(m1+m2)
                  if (j3 < ABS(m3)) cycle
                  jm(1,1)=j1; jm(2,1)=j2; jm(3,1)=j3
                  jm(1,2)=m1; jm(2,2)=m2; jm(3,2)=m3
                  ltemp = Reorder3j(jm, rjm)
                  pos = get3jPos1(rjm(1,1),rjm(2,1),rjm(3,1),rjm(3,2))

                  if (pos>jSize) then 
                      print *, 'Outof Index:', pos, rjm(1:3,1),rjm(3,2)
                      cycle
                  end if
                  if (jLen(pos) == 0)  &
                      jLen(pos)=min(rjm(1,1)-rjm(3,2), rjm(2,1))+rjm(2,1)+1
               end do  ! m2
            end do     ! m1
         end do  ! j3
      end do     ! j2
   end do        ! j1
end

!*****************************************************************************
subroutine get3jIndex0(jmax, jSize, jInd)
   implicit none
   integer, intent(IN) :: jmax(3), jSize
   integer, intent(OUT):: jInd(4, jSize)

   integer :: j1, j2, j3, m1, m2, m3, pos, get3jPos1
   integer :: j3max, j3min, jsum, jm(3,2), rjm(3,2)
   logical :: fill(jSize), ltemp, Reorder3j

   jInd(1:4,1:jSize)=0
   fill(1:jSize) = .FALSE.
   do j1=0, jmax(1)
      do j2 = 0, jmax(2)
         j3min = ABS(j1-j2); j3max = min(jmax(3), j1+j2)
         do j3 = j3min, j3max
            do m1 = -j1, j1
               do m2 = -j2, j2
                  m3 = -(m1+m2)
                  if (j3 < ABS(m3)) cycle
                  jm(1,1)=j1; jm(2,1)=j2; jm(3,1)=j3
                  jm(1,2)=m1; jm(2,2)=m2; jm(3,2)=m3
                  ltemp = Reorder3j(jm, rjm)
                  pos = get3jPos1(rjm(1,1),rjm(2,1),rjm(3,1),rjm(3,2))
                  if ((pos>jSize) .OR. (pos<1)) then
                       print *, 'Error at pos:', pos, rjm
                       cycle
                  end if
                  if (.NOT. fill(pos)) then
                     jInd(1:3,pos) = rjm(1:3,1); jInd(4,pos)=rjm(3,2)
                     fill(pos) = .TRUE.
                  end if
                  
               end do ! m2
            end do    ! m1
         end do     ! j3
      end do        ! j2
   end do           ! j1
end

!*************************************************************************
subroutine get3jLenInd0(jmax, jSize, jLen, jInd)
   implicit none
   integer, intent(IN) :: jmax(3), jSize
   integer, intent(OUT):: jLen(jSize), jInd(4, jSize)

   integer :: j1, j2, j3, m1, m2, m3, pos, get3jPos1
   integer :: j3max, j3min, jsum, jm(3,2), rjm(3,2)
   logical :: fill(jSize), ltemp, Reorder3j

   fill(1:jSize)=.FALSE.
   jLen(1:jSize)=0;    jInd(1:4,1:jSize)=0
   do j1=0, jmax(1)
      do j2 = 0, jmax(2)
         j3min = ABS(j1-j2); j3max = min(jmax(3), j1+j2)
         do j3 = j3min, j3max
            do m1 = -j1, j1
               do m2 = -j2, j2
                  m3 = -(m1+m2)
                  if (j3 < ABS(m3))  cycle
                  jm(1,1)=j1; jm(2,1)=j2; jm(3,1)=j3
                  jm(1,2)=m1; jm(2,2)=m2; jm(3,2)=m3
                  ltemp = Reorder3j(jm, rjm)
                  pos = get3jPos1(rjm(1,1),rjm(2,1),rjm(3,1),rjm(3,2))
                   
                  if (pos>jSize) cycle

                  if (.NOT. fill(pos)) then
                      jInd(1:3,pos) = rjm(1:3,1); Jind(4,pos)=rjm(3,2)
                      jLen(pos)=min(rjm(1,1)-rjm(3,2), rjm(2,1))+rjm(2,1)+1
                      fill(pos) = .TRUE.
                  end if
               end do
            end do
         end do
      end do
   end do
end

!*******************************************************************************
