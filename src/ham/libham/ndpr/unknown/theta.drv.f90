!
! Program for NDP-DVR
!
program test_ndp
   use ndpr
   implicit none
   integer :: nBasis, i, j1,m1,k1,j2,m2,k2, opt
   double precision :: sTh, cTh, sPhi, cPhi
   logical :: callr=.FALSE.
   external :: fitVlr, fitVBR
   double precision :: ct0, ct1

   call CPU_TIME(ct0)
   print *
   print *, ' *****************************************'
   print *, '         Non-Direct-Product DVR '  
   print *, '******************************************'
   print *, ' Read input parameters from STDIN'
   call readNDPRFile('ndp.in')

   opt = 1   
   print *, ' Initializing .........'
   if (init()) then
        if (useSP) then
            print *, ' Reading 1D potential from file:', inFile
            call readVR(inFile)
        else
            if (callr) then
               call initVR(fitVlr)
            else
               call initVR(fitVBR)
            end if
        end if

        do 
           print *, 'Input j1, m1, k1, j2, m2, k2: (j1,j2=[0, n], m1, m2=[0, j1/j2], k1,k2=[-1,1])'
           read(*,*) j1, m1, k1, j2,m2,k2
           call getTheta(j1,m1,j2,m2,sTh, cTh)

           print *
           print *, ' <cosTheta>:', cTH
           cTH = getThetaElem(j1,m1,j2,m2)
           print *, ' <cosPhi>:', cTh

           print *, 'Input option to continue,0 to exit'
           read (*, *) opt
           if (opt==0) exit
        end do
    
   else
        print *, ' Error in initializing.'
   end if 
   call CPU_TIME(ct1)

   print *, ' ======== Finished Program: CPU Time:', ct1-ct0, '======='
   print *
   call final()

10 Format(4F15.9)
 
end
