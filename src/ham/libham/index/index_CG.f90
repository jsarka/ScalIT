!
! Index for Clebsch-Gordon coefficients
!
integer function getCGSize(j1max,j2max,jmax,kmax) 
   implicit none
   integer, intent(IN)  :: j1max, j2max, jmax, kmax 
      
   integer :: jNum, i, j1,j2,j
   integer :: getCGSize1, getCGjSize
   integer, allocatable :: jSize(:)

   jNum = getCGjSize(j1max, j2max, jmax)

   allocate(jSize(jNum))    
   getCGSize = getCGSize1(j1max,j2max,jmax,kmax,jNum,jSize)
   deallocate(jSize)

end

!*********************************************************************
integer function getCGSize1(j1max,j2max,jmax,kmax,N,jSize)
   implicit none
   integer, intent(IN)  :: j1max, j2max, jmax, kmax, N
   integer, intent(OUT) :: jSize(N)
      
   integer :: jNum, i, j, j1, j2
   integer :: getCGjSize, getCGmkSize
   integer,allocatable :: jInd(:,:),jLen(:)

   jNum = getCGjSize(j1max, j2max, jmax);
   
   allocate(jInd(3,jNum), jLen(jNum))
   
   call getCGjIndex(j1max,j2max,jmax,jNum,jInd)  

   DO i=1,jNum
      jLen(i) = getCGmkSize(kmax,jInd(1,i),jInd(2,i),jInd(3,i))
   END DO   

   getCGSize1 = SUM(jLen(1:jNum))
   i = min(N,jNum)
   jSize(1:i) = jLen(1:i)  
   deallocate(jInd, jLen)
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getCGIndex(j1max,j2max,jmax,kmax,N,jmkInd)
   implicit none
   integer, intent(IN) :: j1max,j2max,jmax,kmax,N
   integer, intent(OUT):: jmkInd(5,N)   !(j1j2j,m1(K-m1)K)

   integer :: j1,j2,j, m1, k
   integer :: i0, k0, m2, j0min,j0max

   i0 = 1
   do j1 = 0, j1max
      do j2 = 0, j2max
         j0min = ABS(j1-j2); j0max=min(j1+j2, jmax)
         do j = j0min, j0max
            k0 = min(j, kmax)
            do k = 0, k0
               do m1 = -j1, j1
                  m2 = ABS(k-m1)
                  if (.NOT. (m2>j2)) then
                      jmkInd(1,i0)=j1; jmkInd(2,i0)=j2
                      jmkInd(3,i0)=j;  jmkInd(4,i0)=m1
                      jmkInd(5,i0)=k;  i0 = i0 + 1
                      if (i0 > N) return
                  end if
               end do
            end do
         end do
      end do
   end do
end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function getCGPos(j1max,j2max,jmax,kmax, jmk)
   implicit none
   integer, intent(IN) :: j1max,j2max,jmax,kmax, jmk(5)   

   integer :: j1,j2,j, m1, k
   integer :: i0, k0, m2,j0min,j0max

   i0 = 1
   do j1 = 0, j1max
      do j2 = 0, j2max
         j0min = ABS(j1-j2); j0max=min(j1+j2, jmax)
         do j = j0min, j0max
            k0 = min(j, kmax)
            do k = 0, k0
               do m1 = -j1, j1
                  m2 = ABS(k-m1)
                  if (.NOT. (m2>j2)) then
                     if ((jmk(1)==j1) .AND. (jmk(2)==j2) .AND. &
                         (jmk(3)==j)  .AND. (jmk(4)==m1) .AND. &
                         (jmk(5)==k)) then
                         getCGPos = i0; return
                     else
                        i0 = i0 + 1
                     end if             
                  end if
               end do
            end do
         end do
      end do
   end do
end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getCGInd(j1max,j2max,jmax,kmax,Ind, jmKInd)
   implicit none
   integer, intent(IN) :: j1max,j2max,jmax,kmax,ind
   integer, intent(OUT):: jmkInd(5)

   integer :: j1,j2,j, m1, k
   integer :: i0, k0, m2, j0min,j0max

   i0 = 1
   do j1 = 0, j1max
      do j2 = 0, j2max
         j0min = ABS(j1-j2); j0max=min(j1+j2, jmax)
         do j = j0min, j0max
            k0 = min(j, kmax)
            do k = 0, k0
               do m1 = -j1, j1
                  m2 = ABS(k-m1)
                  if (.NOT. (m2>j2)) then
                     if (i0==ind) then
                        jmkInd(1)=j1; jmkInd(2)=j2
                        jmkInd(3)=j;  jmkInd(4)=m1
                        jmkInd(5)=k;  return
                     else
                        i0 = i0 + 1
                     end if             
                  end if
               end do
            end do
         end do
      end do
   end do
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



