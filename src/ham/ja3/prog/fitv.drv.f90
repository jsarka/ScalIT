!
!  Testing the fitting function for 1D potential
!
program test_fit
   implicit none

   double precision :: r(1), vBr(1), vlr(1)
   integer :: opt = 1

   print *, '*********************************************'
   print *, '  Testing fitting function for 1D potential  '
   print *, '*********************************************' 
  
   do while (opt /= 0) 
      print *, ' Input r value:'
      read *, r(1)
      call fitVBR(1, r, VBR)
      call fitVlr(1, r, vlr)
      write (*, 100) r, vbr, vlr
      print *
      print *, '  Input option to continue: 1: exit / other: continue'
      read *, opt
   end do
  100 format('r=',F10.6, 2x, 'VBR=',F10.6,2x,'vlr=',F10.6)
   print *, ' **************   Finish the program   **************'
end 

