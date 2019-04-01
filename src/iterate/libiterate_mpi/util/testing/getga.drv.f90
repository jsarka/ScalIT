!cccccccccccccccccccccccccccccccccccccc
!c   test program for get***Index     c
!cccccccccccccccccccccccccccccccccccccc
program test_getmindex
   implicit none
   include 'mpif.h'
   integer, parameter :: NMAX=20
   integer :: sF, sN(NMAX)
   integer :: sNode,id, ierr

   integer(kind=MPI_OFFSET_KIND) :: sPos1,ePos1, nDIm1
   integer :: nLen1, nout, nin, bnum1, mDim1
   integer :: i, j

   double precision, allocatable :: x0(:), y0(:), z0(:)

   call MPI_Init(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, id, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, sNode, ierr)

   if (id ==0 ) then
!      print *, ' Input number of degrees of freedom: N='
      read(*,*) sF
   end if
   call MPI_BCAST(sF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
   if (id==0) then
!      print *, ' Input size for each degree of freedom: sN(1:N):'
      read(*,*) sN(1:sF)
   end if
   call MPI_BCAST(sN,sF,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

   nLen1 = 1
   if (id==0) then
      do i = 1, sF
         nLen1 = nLen1*sN(i)
      end do
   end if
   allocate(y0(nLen1), z0(nLen1))
   do i=1,nLen1
        z0(i)=i
   end do
!   z0(1)=1.5D0

   do i = 1, sF
         call getSimplePos(i,sF,sN,sNode,id,nDim1,sPos1,ePos1,bNum1,mDim1)
         nin = ePos1 - sPos1+1
         if (bNum1==0) then
            nout=sN(i)
         else
            nout = 1 
         end if
         allocate(x0(nin*nout))
         call InitGA(i,sF, sN,sNode,id,Nin,nout,x0)
         call getGA(i,sF,sN,sNode,id,Nin*Nout,x0,nLen1,y0)
         if (id==0) then
            print *
            print *, ' Global data: If there is no erro, it is OK at layer=',i
            do j=1,nlen1
               if (y0(j)/=z0(j)) print *, 'error data for index=', j
            end do
         end if
         deallocate(x0)
   end do

   if (id==0) print *, 'finish!'

   deallocate(y0,z0)

   call MPI_FINALIZE(ierr)

end program
