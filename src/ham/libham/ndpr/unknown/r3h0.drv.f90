!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Combine all H0 from different files into 1 file      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program H0Combining
   implicit none
   integer, parameter :: FILENAME_LEN = 128
   integer, parameter :: DOFMAX = 20
   integer :: dof
   integer :: sN(DOFMAX)    
   double precision, allocatable :: H0(:), H1(:) 
   logical :: readMode, saveMode 
   character(LEN=FILENAME_LEN) H0File1, H0File2, H0File
   integer :: nSize, ind , info
   double precision :: ct1, ct2

   call CPU_Time(ct1)
   dof = 2
   read(*,*) sN(1:dof)
   read(*,*) readMode
   read(*, '(A)') H0File1
   read(*, '(A)') H0File2
   read(*, *) saveMode
   read(*, '(A)') H0File

   print *, ' ======================================================='
   print *, '      Combine H0 into one file for A3 in NDPR'
   print *, ' ======================================================='

   print *, ' # of DVR lr:', sN(1), '  # of DVR Br:', sN(2)
   if (readMode) then
      print *, ' Input data are stored in Binary Format'
   else
      print *, ' Input data are stored in ASCII Format '
   end if
   
   nSize = sN(1)**2+sN(2)**2
   ind=1
   allocate(H1(nSize), stat = info)
   if (info == 0) then
      allocate(H0(sN(1)**2))
      print *, '   Reading H0 Data of lr from :', H0File1
      call loadData(sN(1)**2,H0,readMode,H0File1)
      H1(1:sN(1)**2) = H0(1:sN(1)**2)
      deallocate(H0)
      ind = ind + sN(1)**2

      allocate(H0(sN(2)**2))
      print *, '   Reading H0 Data of Br from :', H0File2
      call loadData(sN(2)**2,H0,readMode,H0File2)
      H1(ind:ind+sN(2)**2-1)=H0(1:sN(2)**2)
      deallocate(H0)

      if (saveMode) then  
          print *, ' Save combined H0 in Binary format in File:',H0File
      else
          print *, ' Save combined H0 in ASCII format in File:', H0File
      end if
 
      call saveData(nSize,H1,saveMode,H0File)
     
      deallocate(H1)
 
   else
       print *, ' Error in Allocating Memory'

   end if

   call CPU_Time(ct2)
   print *, ' Total CPU Time:', ct2-ct1
   print *, ' Finish the program!'
   print *, '=============================================='

end 


