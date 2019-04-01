!
! MPI 2 version of ja4 module
! Need Parallel IO
!

program test_mja4
   use mja4
   implicit none
   double precision :: ct1, ct2, ht1, ht2
   
   if (MInit()) then
      if (myid==rootID) then
         call CPU_Time(ct1)
         print *, '*********************************************'
         print *, '     Calculate H0 for Tetratomic Molecule     '
         print *, '*********************************************' 
         print *, '     Read Input Paramaters from STD Input    '
      
         call printJA4()
      end if

      call CPU_Time(ht1)
      call McalSaveH0()
      call CPU_Time(ht2)
      if (myID==rootID) print *, ' CPU Time for H0:', ht2-ht1

      call CPU_Time(ht1)
      call McalSaveHGM()
      call CPU_Time(ht2)
      if (myID==rootID) print *, ' CPU Time for HGM:', ht2-ht1

      if (myID==rootID) then
         call CPU_Time(ct2)
         print *
         print *, ' CPU Time for the program:', ct2-ct1
         print *, '************     Finish the Program ****************'
         print *
      end if
   else
       print *, '  Error in allocating memory in Init()'
   end if

   call MFinal()

end 

