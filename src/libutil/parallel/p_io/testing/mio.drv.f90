program test_Mfile

   implicit none
   include "mpif.h"
   character(LEN=*),parameter  :: fileName = '/pvfs/scratch/wchen/tmp/tmp1.dat'
   integer, parameter :: NIN = 10
   integer, parameter :: dbSize  = 8
   double precision :: buf1(NIN),buf2(NIN)
   integer, parameter :: rootID=0
   integer :: ierr, id, nproc
   integer(KIND=MPI_OFFSET_KIND) :: pos
   integer :: fh, i

   call MPI_Init(ierr)
  
   call MPI_COMM_Rank(MPI_COMM_WORLD, id, ierr)
   call MPI_Comm_Size(MPI_COMM_WORLD, nproc, ierr)
   if (id==0)  print *, ' Initialize MPI: Open file:', filename

   do i = 1, NIN
      buf1(i)=(id*1000+i)
   end do

   pos = id * NIN + 1

   print *, ' Save data at id=', id,' at position= ', pos

   if (id==rootID) print *, 'Data at root :', buf1

   call MSaveData(MPI_COMM_WORLD,filename,pos,NIN,buf1,ierr)

   call MLoadData(MPI_COMM_WORLD,filename,pos,NIN,buf2,ierr)

   if (id==rootID) print *, ' Read data at root:', buf2, buf2-buf1

  call MPI_Finalize()

end  

