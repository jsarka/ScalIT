program test_Mfile

   implicit none
   include "mpif.h"
   character(LEN=*),parameter  :: fileName = '/home/wchen/tmp/tmp_dist.dat'
   integer, parameter :: NIN = 6, NOUT=4
   integer :: N
!   double complex :: buf1(NIN,NOUT),buf2(NIN,NOUT)
   double precision :: buf1(NIN,NOUT),buf2(NIN,NOUT)
   integer, parameter :: rootID=0
   integer :: ierr, id, nproc
   integer(KIND=MPI_OFFSET_KIND) :: pos, gSize
   integer ::  i, j, k

   call MPI_Init(ierr)
  
   call MPI_COMM_Rank(MPI_COMM_WORLD, id, ierr)
   call MPI_Comm_Size(MPI_COMM_WORLD, nproc, ierr)
   if (id==0)  print *, ' Initialize MPI: Open file:', filename, nproc

   if (id==rootID) print *, 'Data at root :', nproc

   do i = 1, NOUT
!       buf1(i)=dcmplx((id*1000+i), i*1000+id)
      do j = 1, NIN
            buf1(j,i)=(id*10000+i*100+j)
!          buf1(j,i)=dcmplx((id*10000+i*100+j), (j*10000+i*100+id))
      end do
   end do

!   if (id==rootID) print *, 'Data at root :', nproc

   N = nproc*nin

!   if (id==rootID) print *, 'call MGather:', buf1

!   call MGatherSeq_CX(nproc, id, N,filename,nin,buf1)
!   call MGatherGrid_CX(nproc, id, N,filename,nin,nout,buf1)
    call MGatherGrid(nproc, id, N,filename,nin,nout,buf1)

   if (id==rootID) print *, ' Finish Gather'
!   n = nproc*Nin

!   buf2(1:nin)=0.0D0;
!   call MScatterSeq_CX(nproc,id,N,filename,nin,buf2)
!   call MScatterGrid_CX(nproc,id,N,filename,nin,nout,buf2)
  call MScatterGrid(nproc,id,N,filename,nin,nout,buf2)

   if (id==rootID) print *, ' diff data at root:', buf2-buf1

  call MPI_Finalize()

end  

