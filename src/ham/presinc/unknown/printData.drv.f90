!
! print 1D data
!
program printData
   implicit none
   integer :: N1,opt
   double precision, allocatable :: data(:)
   integer :: datLen
   character(len=128) :: fname

   print *, ' Input the size of data:'
   read(*,*) N1
   print *, ' Data storage format: <=0: direct access, >0:sequential access'
   print *, ' Input data storage format:'
   read(*,*) opt 
   print *, " Input the file name:starts/ends with ': "
   read(*,*) fname

   allocate(data(N1))

   inquire(IOLENGTH=datLen) data    

   if (opt<=0) then
      open(99,FILE=fName, STATUS='Old', FORM='UNFORMATTED', &
               ACCESS='direct',RECL=datLen)
      read(99, rec=1) data
      close(99)
   else
      open(99,FILE=fName, STATUS='Old', FORM='UNFORMATTED')
      read(99) data(1:N1)
      close(99)
   end if

   print *, ' Data from binary file:', fname

   print *, data(1:N1)

   deallocate(data)


end program

