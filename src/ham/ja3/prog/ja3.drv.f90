program test_ja3
   use ja3
   implicit none
   double precision :: ct1, ct2 
   
   call CPU_Time(ct1)
   print *, '*********************************************'
   print *, '     Calculate H0 for Triatomic Molecule     '
   print *, '*********************************************' 
   print *, '     Read Input Paramaters from STD Input    '
   
   if (initJA3()) then
      call printJA3()

      if ( calSaveH0()) then
         call calSaveHGM()
      else
         print *, ' Error in calculate H0'
      end if
   else
       print *, '  Error in allocating memory in Init()'
   end if

   call finalJA3()
   call CPU_Time(ct2)
   print *
   print *, ' CPU Time for the program:', ct2-ct1
   print *, '************     Finish the Program ****************'
   print *

end 

