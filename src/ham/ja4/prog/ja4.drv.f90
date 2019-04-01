program test_ja4
   use ja4
   implicit none
   double precision :: ct1, ct2 
   
   call CPU_Time(ct1)
   print *, '*********************************************'
   print *, '     Calculate H0 for Tetratomic Molecule     '
   print *, '*********************************************' 
   print *, '     Read Input Paramaters from STD Input    '
   
   if (initJA4()) then
      call printJA4()

      if( calSaveH0()) then
          call calSaveHGM()
      else
          print *, ' Error in calculating H0:'
      end if
   else
       print *, '  Error in allocating memory in Init()'
   end if

   call finalJA4()
   call CPU_Time(ct2)
   print *
   print *, ' CPU Time for the program:', ct2-ct1
   print *, '************     Finish the Program ****************'
   print *

end 

