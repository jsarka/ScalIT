!
!  Testing the fitting function for 1D potential
!
program test_fit
   implicit none

   double precision :: lr1, lr2, Br, theta1,theta2,phi, tmp
   integer :: opt = 1
   double precision :: pot

   print *, '*********************************************'
   print *, '  Testing Potential function for 4-atom      '
   print *, '      Molecule in Jacobi Coordinate          '
   print *, '*********************************************' 
  
   do while (opt /= 0) 
      print *, ' Input radial values: lr1, lr2, Br '
      read *, lr1, lr2, Br
      print *, ' Input Angle Values: theta1, theta2, phi'
      read *, theta1, theta2, phi
      tmp = pot(lr1, lr2, Br, theta1, theta2, phi)
      write (*, 100) lr1, lr2, Br
      write (*, 110) theta1, theta2,phi, tmp
      print *
      print *, '  Input option to continue: 1: exit / other: continue'
      read *, opt
   end do

  100 format(' lr1=',F10.6, 2x, 'lr2=',F10.6,2x,' Br=',F10.6)
  110 format(' theta1=',F10.6, 2x, 'theta2=',F10.6,2x,'phi=',F10.6,2x,' Pot=',F10.6)
   print *, ' **************   Finish the program   **************'
end 

