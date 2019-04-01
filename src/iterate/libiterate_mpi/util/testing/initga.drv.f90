!
! test program for getPos
!
program test_initga
   include 'mpif.h'
   integer, parameter :: NMAX=20
   integer :: sF, sNode, sN(NMAX)
   
   integer :: i, bnum1, mDim1, nin, nout
   integer(kind=MPI_OFFSET_KIND) :: sPos1, ePos1, nDim1
   double precision, allocatable :: x0(:)

 !  print *, ' Input number of the degrees of freedom:'
   read(*,*) sNode, sF
   read(*,*) sN(1:sF)

   print *
   print *, ' # of Layers:',sF
   print *, ' Layer configuration:'
   print *, sN(1:sF)
   print *, ' # of nodes:', sNode

   do i = 1, sF
      print *
      print *, ' Configure for layer i=', i
      do j=1, sNode
         call getSimplePos(i,sF,sN,sNode,j-1,nDim1,sPos1,ePos1,bNum1,mDim1)
         nin = ePos1 - sPos1+1
         if (bNum1==0) then
            nout=sN(i)
         else
            nout = 1 
         end if
         allocate(x0(nin*nout))
         call InitGA(i,sF, sN,sNode,j-1,Nin,nout,x0)
         print *, ' Data for node=', j-1
         print *, 'nDIm=',nDim1,'. sPos=',sPos1,'. ePos=',ePos1
         print *, 'bNum=',bNum1, '. mDim1=',mDim1, ' Local data:'
         print *, x0
         deallocate(x0)
      end do
       
   end do

   print *, 'finish!'


end program
