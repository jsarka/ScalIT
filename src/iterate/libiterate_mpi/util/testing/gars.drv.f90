!cccccccccccccccccccccccccccccccccccccccccccccc
!c   testing program for direct send/recv     c
!cccccccccccccccccccccccccccccccccccccccccccccc
program test_gasr
   implicit none
   include 'mpif.h'
   integer, parameter :: NMAX=20
   integer :: sF, sN(NMAX)
   integer :: sNode,id, ierr

   integer :: level1, level2
   integer(kind=MPI_OFFSET_KIND) :: sPos1,ePos1,nDim1,sPos2,ePos2,nDim2
   integer :: nLen1, bnum1, mDim1, nLen2, bnum2, mDim2
   integer :: nin1, nout1, nin2, nout2
   integer :: i, j, k

   double precision, allocatable :: x1(:),x2(:),y0(:),z0(:)

   call MPI_Init(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, id, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, sNode, ierr)

   if (id ==0 ) then
      read(*,*) sF
   end if
   call MPI_BCAST(sF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
   if (id==0) then
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
!   y0=0.0D0; z0=0.0

   if (id==0) then
       print *, ' Configure: sF=',sF
       print *, ' Grid Size:', sN(1:sF)
   end if

   do i = 1, sF
      do j = 1, sF
        level1 = i; level2 = j
       call getSimplePos(level1,sF,sN,sNode,id,nDim1,sPos1,ePos1,bNum1,mDim1)
       call getSimplePos(level2,sF,sN,sNode,id,nDim2,sPos2,ePos2,bNum2,mDim2)

       if (id==0) then
         print *
         print *, 'level1=', level1, '  level2=',level2
       end if

       nin1 = ePos1 - sPos1+1
       if (bNum1==0) then
          nout1=sN(level1)
       else
          nout1 = 1 
       end if

       nin2 = ePos2 - sPos2+1
       if (bNum2==0) then
          nout2=sN(level2)
       else
          nout2 = 1 
       end if

       allocate(x1(nin1*nout1), x2(nin2*nout2))
       x2=0.0D0
       call InitGA(level1,sF, sN,sNode,id,Nin1,nout1,x1)
       call directRS(sF,sN,sNode,id,level1,Nin1*Nout1,x1,level2,nin2*nout2,x2)
       call getGA(level2,sF,sN,sNode,id,Nin2*Nout2,x2,nLen1,y0)
       if (id==0) then
          print *
          print *, ' Global data: If there is no erro, it is OK between ',level1, level2
!          print *, y0
          do k=1,nlen1
             if (y0(k)/=z0(k)) then
                print *, 'error data for index=', j, ', K=',k,y0(k),z0(k)
           !     eFlag=.true.
             end if
          end do
          
       end if
       deallocate(x1, x2)
       end do
   end do

   if (id==0) print *, 'finish!'

   deallocate(y0, z0)

   call MPI_FINALIZE(ierr)

end program
