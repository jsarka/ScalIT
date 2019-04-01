!
! Program for NDP-DVR
!
program test_ndp
   use ndpr
   implicit none
   integer :: nBasis, i, mmax, j1,m1,k1,j2,m2,k2, opt
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
           print *, 'Input max m:'
           read(*,*) mmax
           do m1 = 0, mmax
              do m2 = 0, mmax
                 k1 = 1;  k2 = 1
                 call getFmk(m1,k1,m2, k2, sPhi, cPhi)
                 write(*,10) m1,k1,m2, k2, sPhi, cPhi

                 k1 = 1;  k2 = -1
                 call getFmk(m1,k1,m2, k2, sPhi, cPhi)
                 write(*,10) m1,k1,m2, k2, sPhi, cPhi

                 k1 = -1;  k2 = 1
                 call getFmk(m1,k1,m2, k2, sPhi, cPhi)
                 write(*,10) m1,k1,m2, k2, sPhi, cPhi

                 k1 = -1;  k2 = -1
                 call getFmk(m1,k1,m2, k2, sPhi, cPhi)
                 write(*,10) m1,k1,m2, k2, sPhi, cPhi

              end do
           end do

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

10 Format('m1=',I4,1x,'k1=',I2,1x,'m2=',I4,1x,'k2=',I2,1x,'sPhi:',F10.7,1x,'cPhi:',F10.7)
 
end
