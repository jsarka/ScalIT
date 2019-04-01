program test_Mfile

   implicit none
   include "mpif.h"
   character(LEN=*),parameter  :: fileName = '/pvfs/scratch/wchen/tmp/tmp3.dat'
   integer, parameter :: NIN = 3,NOUT=2
   integer, parameter :: dbSize  = 8
   double precision :: buf1(NIN,NOUT,NOUT),buf2(NIN,NOUT,NOUT)
   integer, parameter :: rootID=0
   integer :: ierr, id, nproc
   integer(KIND=MPI_OFFSET_KIND) :: pos, gSize
   integer ::  i, j, k

   call MPI_Init(ierr)
  
   call MPI_COMM_Rank(MPI_COMM_WORLD, id, ierr)
   call MPI_Comm_Size(MPI_COMM_WORLD, nproc, ierr)
   if (id==0)  print *, ' Initialize MPI: Open file:', filename

   do i = 1, NOUT
      do j = 1, NOUT
         do k = 1, NIN
            buf1(k,j,i)=(id*10000+i*1000+j*100+k)
         end do
      end do
   end do

   pos = id * NIN + 1; gSize=nproc*NIN

   print *, ' Save data at id=', id,' at position= ', pos

   if (id==rootID) print *, 'Data at root :', buf1

   call MSaveDataGrid(MPI_COMM_WORLD,filename,pos,gSize,NIN,NOUT,buf1,ierr)

   call MLoadDataGrid(MPI_COMM_WORLD,filename,pos,gSize,NIN,NOUT,buf2,ierr)

   if (id==rootID) print *, ' Read data at root:', buf2, buf2-buf1

  call MPI_Finalize()

end  

