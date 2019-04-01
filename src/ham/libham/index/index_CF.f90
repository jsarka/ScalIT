!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Handle indices for centrifugal potential problem    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function getCFSize(jmax, jSize)
   implicit none
   integer, intent(IN) :: jmax, jSize(jmax+1)

   integer :: j

   getCFSize = 0
   do j=0, jmax
      getCFSize = getCFSize + (j+j+1)*jSize(j+1)
   end do

end


!cccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getCFIndex(jmax, jSize, N, ind)
   implicit none
   integer, intent(IN) :: jmax, jSize(jmax+1),N
   integer, intent(OUT):: ind(4,N)

   integer :: j, jj, m, i0

   i0 = 1
   do j = 0, jmax
      do jj = 1, jSize(j+1)
         ind(1, i0) = j; ind(2, i0) = jj
         ind(3, i0) = 0; ind(4, i0) = 1
         i0 = i0 + 1
         if (i0 > N) return
         do m = 1, j
            ind(1, i0) = j; ind(2, i0) = jj
            ind(3, i0) = m; ind(4, i0) = 1 
            i0 = i0 + 1
            if (i0 > N) return
      
            ind(1, i0) = j;  ind(2, i0) = jj
            ind(3, i0) = m;  ind(4, i0) = -1 
            i0 = i0 + 1
            if (i0 > N) return
         end do
      end do
   end do
end
!ccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function getCFPos(jmax, jSize, ind)
   implicit none
   integer, intent(IN) :: jmax, jSize(jmax+1), ind(4)

   integer :: j, jj, m, i0 

   getCFPos = 0; i0 = 1
   do j = 0, jmax     
      do jj = 1, jSize(j+1)
         if ((ind(1)==j).AND.(ind(2)==jj).AND.(ind(3)==0)) then
             getCFPos = i0; return
         else
             i0 = i0 + 1 
         end if
         do m = 1, j
            if ((ind(1)==j).AND.(ind(2)==jj).AND.(ind(3)==m)) then
                if (ind(4)==1) then
                   getCFPos = i0  
                else
                   getCFPos = i0 + 1
                end if       
                return
            else
               i0 = i0 + 2
            end if          
         end do
      end do
   end do
end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getCFInd(jmax, jSize, pos, ind)
   implicit none
   integer, intent(IN) :: jmax, jSize(jmax+1),pos
   integer, intent(OUT):: ind(4)

   integer :: j, jj, m, i0

   i0 = 1
   do j = 0, jmax
      do jj = 1, jSize(j+1)
         if (pos==i0) then
             ind(1)=j; ind(2)=jj; ind(3)=0;ind(4)=1; return
         else
             i0 = i0 + 1
         end if
         do m = 1, j
            if (pos==i0) then
                ind(1)=j; ind(2)=jj; ind(3)=m;ind(4)=1; return
            else
                if ((i0+1)==pos) then
                   ind(1)=j;ind(2)=jj;ind(3)=m;ind(4)=-1;return
                else
                   i0 = i0 + 2
                end if
            end if
         end do
      end do
   end do
end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc






