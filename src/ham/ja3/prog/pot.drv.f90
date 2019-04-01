!
!  Testing the fitting function for 1D potential
!
program test_pot
   implicit none
   character(LEN=*),parameter  :: VBRFILE = 'vBR.dat'
   character(LEN=*), parameter :: vlrFile = 'vlr.dat'
   double precision, parameter :: lrmin=0.8D0, lrmax=3.5D0
   double precision, parameter :: BRMin=0.0001D0, BRMax=3.5D0
   integer, parameter :: N = 3001
   logical, parameter :: saveMode = .false.

   double precision :: r0(N), v0(N), dr
   integer :: i

   print *, '*********************************************'
   print *, '  Testing fitting function for 1D potential  '
   print *, '*********************************************' 
   print *, ' BRmin:', BRMin, '   BRMax:',BRMAX, ' Num:', N
  
   dr = (BRMax-BRMin)/(N-1)
   do i = 1, N
      r0(i)= (i-1)*dr+BRMin
   end do
   print *, ' Calculate VBR Potential using Fitting Function'
   call fitVBR(N, r0, V0)
   print *, '  Save calculated potential in ', VBRFILE
   call save2Data(N, r0, V0, saveMode, VBRFILE)

   print *, ' lrmin:', lrmin, '   lrmax:',lrmax, ' Num:', N
   dr = (lrMax-lrMin)/(N-1)
   do i = 1, N
      r0(i)= (i-1)*dr+lrMin
   end do
   print *, ' Calculate Vlr Potential using Fitting Function'
   call fitVlr(N, r0, V0)
   print *, '  Save calculated potential in ', vlrfile
   call save2Data(N, r0, V0, saveMode, vlrfile)
   print *, ' **************   Finish the program   **************'
end 

