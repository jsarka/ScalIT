!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                  Distribute the data                    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MDistIndex(nNodes, N1, pInd, pLen)
   implicit none
   integer, intent(IN)  :: nNodes, N1
   integer, intent(OUT) :: pInd(nNodes), pLen(nNodes)

   integer :: pnum, qnum, i

   pNum = N1 / nNodes
   qNum = N1 - pNum * nNodes

   do i = 0, nNodes-1
      if (i<qNum) then
         pLen(i+1) = pNum + 1
      else
         pLen(i+1) = pNum
      end if
   end do 

   pInd(1)=1
   do i = 1, nNodes-1
      pInd(i+1) = pInd(i)+pLen(i)
   end do

end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccc    The Root Distribute the data  ccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MScatterSeq( nNodes, id, N1, fname, N2, data1)
   implicit none
   include 'mpif.h'
   integer, intent(IN)  :: nNodes, id, N1, N2
   character(len=*), intent(IN)  :: fname
   double precision, intent(OUT) :: data1(N2)

   integer, parameter :: DAT_TAG = 10
   double precision, allocatable :: gData(:)
   integer :: pInd(nNodes), pLen(nNodes), ierr, i
   integer :: status(MPI_STATUS_SIZE)


   call MDistIndex(nNodes, N1, pInd, pLen)

   if (id==0) then

       allocate(gData(N1))

       call loadDataDir(N1, gdata, fname)

       data1(1:pLen(1)) = gdata(1:pLen(1))

       do i = 1, nNodes-1
          call MPI_Send(gData(pInd(i+1)),pLen(i+1),MPI_DOUBLE_PRECISION, &
                        i, DAT_TAG, MPI_COMM_WORLD, ierr)          
       end do 

       deallocate(gData)

   else

       call MPI_Recv(data1,pLen(id+1),MPI_DOUBLE_PRECISION, &
                      0, DAT_TAG, MPI_COMM_WORLD, status, ierr)

   end if

end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MScatterSeq_CX( nNodes, id, N1, fname, N2, data1)
   implicit none
   include 'mpif.h'
   integer, intent(IN)  :: nNodes, id, N1, N2
   character(len=*), intent(IN)  :: fname
   double complex, intent(OUT) :: data1(N2)

   integer, parameter :: DAT_TAG = 20
   double complex, allocatable :: gData(:)
   integer :: pInd(nNodes), pLen(nNodes), ierr , i
   integer :: status(MPI_STATUS_SIZE)


   call MDistIndex(nNodes, N1, pInd, pLen)

   if (id==0) then

       allocate(gData(N1))
      
       call loadDataDir_CX(N1, gdata, fname)

       data1(1:pLen(1)) = gdata(1:pLen(1))

       do i = 1, nNodes-1
          call MPI_Send(gData(pInd(i+1)),pLen(i+1),MPI_DOUBLE_COMPLEX, &
                        i, DAT_TAG, MPI_COMM_WORLD, ierr)          
       end do 

       deallocate(gData)
       
   else

       call MPI_Recv(data1,pLen(id+1),MPI_DOUBLE_COMPLEX, &
                      0, DAT_TAG, MPI_COMM_WORLD, status, ierr)

   end if

end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MScatterGrid( nNodes, id, N1, fname, N2, N3, data1)
   implicit none
   include 'mpif.h'
   integer, intent(IN)  :: nNodes, id, N1, N2, N3
   character(len=*), intent(IN)  :: fname
   double precision, intent(OUT) :: data1(N2, N3)

   integer, parameter :: DAT_TAG = 30
   double precision, allocatable :: gData(:,:)
   integer :: pInd(nNodes), pLen(nNodes), ierr   
   integer :: status(MPI_STATUS_SIZE)
   integer :: i, j

   call MDistIndex(nNodes, N1, pInd, pLen)


   if (id==0) then

       allocate(gData(N1, N3))

       call loadDataDir(N1*N3, gdata, fname)

       do j = 1, N3
          data1(1:pLen(1), j) = gdata(1:pLen(1), j)

          do i = 1, nNodes-1
             call MPI_Send(gData(pInd(i+1), j),pLen(i+1),MPI_DOUBLE_PRECISION, &
                        i, DAT_TAG, MPI_COMM_WORLD, ierr)          
          end do 
       end do

       deallocate(gData)

   else

       do j = 1, N3
          call MPI_Recv(data1(1,j),pLen(id+1),MPI_DOUBLE_PRECISION,   &
                      0, DAT_TAG, MPI_COMM_WORLD,status, ierr)
       end do

   end if

end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MScatterGrid_CX( nNodes, id, N1, fname, N2, N3, data1)
   implicit none
   include 'mpif.h'
   integer, intent(IN)  :: nNodes, id, N1, N2, N3
   character(len=*), intent(IN)  :: fname
   double complex, intent(OUT) :: data1(N2, N3)

   integer, parameter :: DAT_TAG = 40
   double complex, allocatable :: gData(:,:)
   integer :: pInd(nNodes), pLen(nNodes), ierr   
   integer :: status(MPI_STATUS_SIZE)
   integer :: i, j

   call MDistIndex(nNodes, N1, pInd, pLen)

   if (id==0) then

       allocate(gData(N1, N3))

       call loadDataDir_CX(N1*N3, gdata, fname)

       do j = 1, N3
          data1(1:pLen(1), j) = gdata(1:pLen(1), j)

          do i = 1, nNodes-1
             call MPI_Send(gData(pInd(i+1), j),pLen(i+1),MPI_DOUBLE_COMPLEX, &
                        i, DAT_TAG, MPI_COMM_WORLD, ierr)          
          end do 
       end do

       deallocate(gData)

   else

       do j = 1, N3
          call MPI_Recv(data1(1,j),pLen(id+1),MPI_DOUBLE_COMPLEX, &
                      0, DAT_TAG, MPI_COMM_WORLD, status, ierr)
       end do

   end if

end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



!ccccccccccc    The root Gather the data from each node  ccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MGatherSeq( nNodes, id, N1, fname, N2, data1)
   implicit none
   include 'mpif.h'
   integer, intent(IN)  :: nNodes, id, N1, N2
   character(len=*), intent(IN) :: fname
   double precision, intent(IN) :: data1(N2)

   integer, parameter :: DAT_TAG = 100
   double precision, allocatable :: gData(:)
   integer :: pInd(nNodes), pLen(nNodes), ierr, i
   integer :: status(MPI_STATUS_SIZE)

   call MDistIndex(nNodes, N1, pInd, pLen)

   if (id==0) then

       allocate(gData(N1))

       gdata(1:pLen(1)) = data1(1:pLen(1))

       do i = 1, nNodes-1
          call MPI_Recv(gData(pInd(i+1)),pLen(i+1),MPI_DOUBLE_PRECISION, &
                        i, DAT_TAG, MPI_COMM_WORLD, status,ierr)          
       end do 

       call saveDataDir(N1, gdata, fname)

       deallocate(gData)

   else

       call MPI_Send(data1,pLen(id+1),MPI_DOUBLE_PRECISION, &
                      0, DAT_TAG, MPI_COMM_WORLD, ierr)

   end if

end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MGatherSeq_CX( nNodes, id, N1, fname, N2, data1)
   implicit none
   include 'mpif.h'
   integer, intent(IN)  :: nNodes, id, N1, N2
   character(len=*), intent(IN)  :: fname
   double complex, intent(IN) :: data1(N2)

   integer, parameter :: DAT_TAG = 200
   double complex, allocatable :: gData(:)
   integer :: pInd(nNodes), pLen(nNodes), ierr , i
   integer :: status(MPI_STATUS_SIZE)


   call MDistIndex(nNodes, N1, pInd, pLen)

   if (id==0) then

       allocate(gData(N1))

       gdata(1:pLen(1)) = data1(1:pLen(1))

       do i = 1, nNodes-1
          call MPI_Recv(gData(pInd(i+1)),pLen(i+1),MPI_DOUBLE_COMPLEX, &
                        i, DAT_TAG, MPI_COMM_WORLD, status, ierr)          
       end do 

       call saveDataDir_CX(N1, data1, fname)

       deallocate(gData)

   else

       call MPI_Send(data1,pLen(id+1),MPI_DOUBLE_COMPLEX, &
                      0, DAT_TAG, MPI_COMM_WORLD, ierr)

   end if

end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MGatherGrid( nNodes, id, N1, fname, N2, N3, data1)
   implicit none
   include 'mpif.h'
   integer, intent(IN)  :: nNodes, id, N1, N2, N3
   character(len=*), intent(IN)  :: fname
   double precision, intent(OUT) :: data1(N2, N3)

   integer, parameter :: DAT_TAG = 300
   double precision, allocatable :: gData(:,:)
   integer :: pInd(nNodes), pLen(nNodes), ierr   
   integer :: i, j
   integer :: status(MPI_STATUS_SIZE)


   call MDistIndex(nNodes, N1, pInd, pLen)

   if (id==0) then

       allocate(gData(N1, N3))

       do j = 1, N3
          gdata(1:pLen(1), j) = data1(1:pLen(1), j)

          do i = 1, nNodes-1
             call MPI_Recv(gData(pInd(i+1), j),pLen(i+1),MPI_DOUBLE_PRECISION, &
                        i, DAT_TAG, MPI_COMM_WORLD,status, ierr)          
          end do 
       end do

       call saveDataDir(N1*N3, gdata, fname)

       deallocate(gData)

   else

       do j = 1, N3
          call MPI_Send(data1(1,j),pLen(id+1),MPI_DOUBLE_PRECISION,   &
                      0, DAT_TAG, MPI_COMM_WORLD, ierr)
       end do

   end if

end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MGatherGrid_CX( nNodes, id, N1, fname, N2, N3, data1)
   implicit none
   include 'mpif.h'
   integer, intent(IN)  :: nNodes, id, N1, N2, N3
   character(len=*), intent(IN)  :: fname
   double complex, intent(OUT) :: data1(N2, N3)

   integer, parameter :: DAT_TAG = 400
   double complex, allocatable :: gData(:,:)
   integer :: pInd(nNodes), pLen(nNodes), ierr   
   integer :: status(MPI_STATUS_SIZE)
   integer :: i, j

   call MDistIndex(nNodes, N1, pInd, pLen)

   if (id==0) then

       allocate(gData(N1, N3))

       do j = 1, N3
          gdata(1:pLen(1), j) = data1(1:pLen(1), j)

          do i = 1, nNodes-1
             call MPI_Recv(gData(pInd(i+1), j),pLen(i+1),MPI_DOUBLE_COMPLEX, &
                        i, DAT_TAG, MPI_COMM_WORLD, status,ierr)          
          end do 
       end do

       call saveDataDir_CX(N1*N3, gdata, fname)

       deallocate(gData)

   else

       do j = 1, N3
          call MPI_Send(data1(1,j),pLen(id+1),MPI_DOUBLE_COMPLEX, &
                      0, DAT_TAG, MPI_COMM_WORLD, ierr)
       end do

   end if

end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
