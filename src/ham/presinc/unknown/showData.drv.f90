!
!
!
program showData
   implicit none
   integer :: N1, M1
   double precision, allocatable :: data(:)
   integer :: datLen
   character(len=128) :: fname

   read(*,*) N1, M1
   read(*,*) fname

   allocate(data(N1))

   inquire(IOLENGTH=datLen) data    

   open(99,FILE=fName, STATUS='Old', FORM='UNFORMATTED', &
               ACCESS='direct',RECL=datLen)
   read(99, rec=1) data
   close(99)

   call sort(N1,data)

   print *, ' Data from binary file:', fname

   print *, data(1:M1)

   deallocate(data)


end program

subroutine sort(N, data)
   implicit none
   integer, intent(IN) :: N
   double precision, intent(INOUT) :: data(N)

   integer :: i, j
   double precision :: tmp

   do i = 1, N
      do j = i+1, N
          if (data(i)>data(j)) then
             tmp=data(i); data(i)=data(j);data(j)=tmp
          end if
      end do
   end do     

end subroutine
