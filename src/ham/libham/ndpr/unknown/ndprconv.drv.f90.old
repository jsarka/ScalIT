!
! Program for NDP-DVR
!
program test_ndpconv
   use ndprconv
   implicit none

   double precision :: ct0, ct1
   integer :: fnum

   call CPU_TIME(ct0)
   print *
   print *, ' *********************************************'
   print *, '  Convergence Test for Non-Direct-Product DVR '  
   print *, '**********************************************'
   print *, ' Read input parameters from STDIN'
   call readInput()
   
   print *, ' Initializing .........'
   if (readVR(VR_FileName)) then
        call printParam()
        print *, ' Reading 1D potential from file:', VR_FileName
        fnum = doConv()
        print *, ' Final N number:', fnum
    else 
        Print *, ' Error in initializing.'
   end if 
   call CPU_TIME(ct1)

   print *, ' ======== Finished Program: CPU Time:', ct1-ct0, '======='
   print *
   call final()

end
