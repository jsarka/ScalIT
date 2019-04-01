program test_NSeq

   implicit none
   include "mpif.h"
   integer, parameter :: FILELEN=128
   character(LEN=FILELEN)  :: fileName 
   integer :: NIN
   double precision,allocatable :: buf1(:),buf2(:)
   integer, parameter :: rootID=0
   integer :: ierr, id, nproc, gSize
   integer(KIND=MPI_OFFSET_KIND) :: pos
   integer ::  i, j

   call MPI_Init(ierr)
  
   call MPI_COMM_Rank(MPI_COMM_WORLD, id, ierr)
   call MPI_Comm_Size(MPI_COMM_WORLD, nproc, ierr)

   if (id==rootID) then
      read(*,'(A)') filename
      read(*,*) NIN
      print *, ' Initialize MPI: Open file:', filename
   end if

   call MPI_BCAST(filename,FILELEN,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(NIN,1,MPI_INTEGER,rootID,MPI_COMM_WORLD,ierr)

   allocate(buf1(NIN),buf2(NIN))

   do j = 1, NIN
      buf1(j)=(id*1000+j)
   end do

   pos = id * NIN + 1; gSize=nproc*NIN

   call NSaveDataSeq(MPI_COMM_WORLD,filename,pos,NIN,buf1,ierr)

   call NLoadDataSeq(MPI_COMM_WORLD,filename,pos,NIN,buf2,ierr)

   print *, ' Data differenceat node=',id,':',  buf2-buf1

   deallocate(buf1)
   deallocate(buf2)

   call MPI_Finalize()

end  

