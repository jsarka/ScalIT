!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Handle indices for centrifugal potential problem    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function getCFjSize(jmax, jSize)
   implicit none
   integer, intent(IN) :: jmax, jSize(jmax+1)

   getCFjSize = SUM(jSize(1:(jmax+1)))
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getCFjIndex(jmax, jSize, N, ind)
   implicit none
   integer, intent(IN) :: jmax, jSize(jmax+1), N
   integer, intent(OUT):: ind(2,N)

   integer :: j, jj, i0

   i0 = 1
   do j = 0, jmax
      do jj = 1, jSize(j+1)
         ind(1, i0) = j; ind(2, i0) = jj         
         i0 = i0 + 1
         if (i0 > N) return
      end do
   end do
end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function getCFjPos(jmax, jSize, ind)
   implicit none
   integer, intent(IN) :: jmax, jSize(jmax+1), ind(2)

   integer :: j, jj, i0

   getCFjPos = 0; i0 = 1
   do j = 0, jmax    
      do jj = 1, jSize(j+1)
         if ((ind(1)==j) .AND. (ind(2)==jj) )then
             getCFjPos = i0
             return
         else
             i0 = i0 + 1
         end if         
      end do
   end do
end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getCFjInd(jmax, jSize, pos, ind)
   implicit none
   integer, intent(IN) :: jmax, jSize(jmax+1),pos
   integer, intent(OUT):: ind(2)

   integer :: j, jj, i0

   i0 = 1
   do j = 0, jmax
      do jj = 1, jSize(j+1)
         if (pos==i0) then
             ind(1)=j; ind(2)=jj;  return
         else
             i0 = i0 + 1
         end if
      end do
   end do
end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




