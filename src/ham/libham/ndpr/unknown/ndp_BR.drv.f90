!
! Program for NDP-DVR
!
program test_ndp_Br
   use ndpr
   implicit none
   integer :: nBasis, i
   double precision, allocatable :: xyz0(:,:), hxyz0(:,:), r0(:)
   external :: fitVBR
   double precision :: ct0, ct1

   call CPU_TIME(ct0)
   print *
   print *, ' *****************************************'
   print *, '         Non-Direct-Product DVR '  
   print *, '******************************************'
   print *, ' Read input parameters from STDIN'
   call readNDPR()
   
   print *, ' Initializing .........'
   if (init()) then

      if (useSP) then
          print *, ' Reading 1D potential from file:', inFile
          call readVR(inFile)
      else
         call initVR(fitVBR)
      end if

      call printParam()

      nBasis = getBasisSize()
      print *, ' Final number of Basis functions:',nBasis
      allocate (xyz0(3, nBasis),r0(nBasis), hxyz0(nBasis, nBasis))
      print *, ' Calculating X and H for DVR ..............'

      if ( getXH_XYZ(xyz0, hxyz0) ) then
         r0(1:nBasis) = DSQRT(xyz0(1,1:nBasis)**2+xyz0(2,1:nBasis)**2+xyz0(3,1:nBasis)**2)
         print *
         print *, ' DVR points: Radius,  (x, y, z)'
         write(*, *) r0(1:nBasis)

         print *
         print *, ' Savin DVR points in file ', xFile
         call saveData(3*nBasis, xyz0, .FALSE., xFile)
 
         print *
         print *, ' Savin HDVR matrix in file ', hFile
         call saveData(nBasis*nBasis, hxyz0, .FALSE., hFile)
      else
         print *, ' Error in calling function getXH_XYZ!'
      end if
      deallocate (xyz0, hxyz0, r0)
   else
      print *, ' Error in initializing.'
   end if 
   call CPU_TIME(ct1)

   print *, ' ======== Finished Program: CPU Time:', ct1-ct0, '======='
   print *
   call final()

10 Format(4F15.9)
 
end
