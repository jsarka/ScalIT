!
! Calculate Residual potential for A3 molecule
!

program NDPRPOT
   implicit none
   integer, parameter :: FILENAME_LEN = 128
   integer :: N1, N2     ! Non-DVR for r and R, r:inner layer, R: outer layer 
   double precision, allocatable :: XYZ1(:,:), XYZ2(:,:), Pot1(:), pot2(:),Res(:,:) 
   logical :: readMode, saveMode 
   character(LEN=FILENAME_LEN) inFile1, inFile2, resFile
   external :: fitVlr, fitVBr, potJA3     ! function to calculate potential
   integer :: i, j, info
   double precision :: ct1, ct2

   call CPU_Time(ct1)
   read(*,*) N1, N2
   read(*,*) readMode
   read(*, '(A)') inFile1
   read(*, '(A)') inFile2
   read(*, *) saveMode
   read(*, '(A)') resFile

   print *, ' ======================================================='
   print *, '      Calculate Residual Potential for A3 in NDPR'
   print *, ' ======================================================='

   print *, ' # of N1:', N1, '  # of N2:', N2
   if (readMode) then
      print *, ' Input data are stored in Binary Format'
   else
      print *, ' Input data are stored in ASCII Format '
   end if
   print *, ' Read N1 Data from ', inFile1
   print *, ' Read N2 data from ', inFile2
   
   allocate(XYZ1(3,N1), XYZ2(3,N2), Pot1(N1), Pot2(N2), Res(N1, N2), stat = info)
   if (info == 0) then
      print *, '   Reading Coordinating Data .......'
      call loadData(3*N1,xyz1,readMode,inFile1)
      call loadData(3*N2,xyz2,readMode,inFile2)

      print *, ' Calculating Potential ........'
      call potFunc3D(N1, XYZ1, fitVlr, pot1)
      print *, 'pot1:', pot1
      print *

      call potFunc3D(N2, XYZ2, fitVBr, pot2)    
      print *, 'pot2:', pot2
      print *

      call potFunc6D(N1, XYZ1, N2, XYZ2, potJA3, res)

      do i = 1, N2
         do j = 1, N1
            RES(j,i) = Res(j,i) - Pot1(j) - Pot2(i)
         end do
      end do
       
      if (saveMode) then
          print *, ' Save RES in Binary format in File:',resFile
      else
          print *, ' Save RES in ASCII format in File:', resFile
      end if
 
      call saveDataDir(N1*N2,Res,resFile)
     
      deallocate(xyz1,xyz2,pot1,pot2,res)
 
   else
       print *, ' Error in Allocating Memory'

   end if

   call CPU_Time(ct2)
   print *, ' Total CPU Time:', ct2-ct1
   print *, ' Finish the program!'
   print *, '=============================================='

end 


