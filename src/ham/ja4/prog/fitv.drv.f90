!
!  Testing the fitting function for 1D potential
!
program test_fit
   implicit none

   double precision :: r(1), vBr(1), vlr1(1), vlr2(1)
   integer :: opt = 1
   external :: fitVlr1, fitVlr2, figVBr

   print *, '*********************************************'
   print *, '  Testing fitting function for 1D potential  '
   print *, '*********************************************' 
  
   do while (opt /= 0) 
      print *, ' Input r value:'
      read *, r(1)
      call fitVBR(1, r, VBR)
      call fitVlr1(1, r, vlr1)
      call fitVlr2(1,r, vlr2)
      write (*, 100) r, vbr, vlr1, vlr2
      print *
      print *, '  Input option to continue: 0: exit / other: continue'
      read *, opt
   end do
  100 format('r=',F10.6, 2x, 'VBR=',F10.6,2x,'vlr1=',F10.6, 2x, 'vlr2=',F10.6)
   print *, ' **************   Finish the program   **************'
end 

