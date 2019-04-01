!
! template program to get wave-function
!
program test_wave3jb
   use wave3jb

   implicit none

   double precision :: ct1, ct2

   print *
   print *, ' *********************************************************'
   print *, ' *   Calculate wave-functions for tri-atomic Molecules   *'
   print *, ' *            Using A-BC Jacobi Coordinator              *'
   print *, ' *********************************************************'

   call cpu_time(ct1)

   if (initW3()) then

      call printW3()

      print *
      print *, ' Calculate & Save Wave Functions ......'
      call calSaveWF()
   else
      print *
      print *, ' Error in initW3A subroutine!'
   end if

   call finalW3()

   call cpu_time(ct2)

   print *
   print *, '  Finish the calculation: CPU Time:', ct2-ct1
   print *, '*******************************************************'
   print *
 
end 


