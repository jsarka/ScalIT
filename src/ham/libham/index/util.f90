!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       Reorder (j1,j2,j3) so (k1>=k2>=k3)           c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine mmm(j1, j2, j3, k1, k2, k3)
   implicit none
   integer, intent(IN)  :: j1, j2, j3
   integer, intent(OUT) :: k1, k2, k3

   k1 = max(j1, j2)
   k3 = min(j1, j2)
   
   if (k1 < j3) then
       k2 = k1; k1=j3
   else
       if ( k3 < j3 ) then
           k2 = j3
       else
           k2 = k3; k3 = j3
       end if
   end if
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function isTri(j1, j2, j3)
   implicit none
   integer, intent(IN)  :: j1, j2, j3
   
   isTri = ((j1>=ABS(j2-j3)) .AND. (j1<=(j2+j3)) .AND. &
            (j2>=ABS(j3-j1)) .AND. (j1<=(j3+j1)) .AND. &
            (j3>=ABS(j1-j2)) .AND. (j1<=(j1+j2)) )
end

logical function isNotTri(j1, j2, j3)
   implicit none
   integer, intent(IN)  :: j1, j2, j3
   
   isNotTri = ((j1<ABS(j2-j3)) .OR. (j1>(j2+j3)) .OR.  &
               (j2<ABS(j3-j1)) .OR. (j1>(j3+j1)) .OR.  &
               (j3<ABS(j1-j2)) .OR. (j1<(j1+j2)) )
end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
