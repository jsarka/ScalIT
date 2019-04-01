!
! Program to test getm3
!

program test_getm3

   integer, dimension(3,6) :: j0 =  &
        RESHAPE( (/ 1,2,3,1,3,2,2,1,3,2,3,1,3,1,2,3,2,1 /), (/3,6/));
   integer :: j1(3), ind(3)
   logical :: getmmm

   DO i=1,6
      print *, 'Original j0:', j0(1:3, i)
      if (getmmm(j0(1,i), j1, ind)) then
         print *, 'Even permutation'
      else
         print *, 'Odd permutation'
      end if

      print *, 'j1:', j1
      print *, 'ind:', ind
   end DO

end
   
   


