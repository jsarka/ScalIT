!
! print 1D data
!
program diffData
   implicit none
   integer :: N1,opt
   double precision, allocatable :: data1(:),data2(:)
   integer :: datLen, i
   double precision :: th, maxErr;
   character(len=128) :: fname1,fname2

   print *, ' Input the size of data:'
   read(*,*) N1
   print *, ' Data storage format: <=0: direct access, >0:sequential access'
   print *, ' Input data storage format:'
   read(*,*) opt 
   print *, ' Input difference threshold:'
   read(*,*) th
   print *, " Input the file name:starts/ends with ': "
   read(*,*) fname1
   read(*,*) fname2

   allocate(data1(N1),data2(N1))

   inquire(IOLENGTH=datLen) data1    

   if (opt<=0) then
      open(99,FILE=fName1, STATUS='Old', FORM='UNFORMATTED', &
               ACCESS='direct',RECL=datLen)
      read(99, rec=1) data1
      close(99)
      open(99,FILE=fName2, STATUS='Old', FORM='UNFORMATTED', &
               ACCESS='direct',RECL=datLen)
      read(99, rec=1) data2
      close(99)
   else
      open(99,FILE=fName1, STATUS='Old', FORM='UNFORMATTED')
      read(99) data1
      close(99)
      open(99,FILE=fName2, STATUS='Old', FORM='UNFORMATTED')
      read(99) data2
      close(99)
   end if

   print *, ' Data from binary files:', fname1,fname2
   data1(1:N1)=ABS(data1(1:N1)-data2(1:N1))
   datLen=0; maxErr=0.0D0
   do i=1, N1
      if (data1(i)>th)  datLen=datLen+1 
      if (data1(i)>maxErr) maxErr=data1(i)
   end do

   print *, ' Data Size:', N1
   print *, ' The difference threshold:', th
   print *, ' Number of difference greater than threshold:', datLen
   print *, ' The maximum difference between the datas:', maxErr

   deallocate(data1,data2)


end program

