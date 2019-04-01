!ccccccccccccccccccccccccccccccccccccccccccccc
!c       IO subroutines in MPI               c
!ccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PSaveDataSeq(myID,rootID,nNodes,locDim,comm,filename,N,data1,ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer, intent(IN)  :: myID,rootID,nNodes
    integer, intent(IN)  :: locDim(nNodes), comm, N
    double precision, intent(IN) :: data1(N)
    integer, intent(OUT) :: ierr
    integer, parameter :: MY_DATA=1000
    
    integer :: dbSize, status(MPI_STATUS_SIZE)
    double precision :: db
    double precision, allocatable :: gOutH(:)
    integer :: disp(nNodes), glen, info 
    integer :: i,j,dstID, recInd
   
    if (myID==rootID) then
        glen=locDim(1)
        do i = 2, nNodes
           if (glen>locDim(i)) glen=locDim(i)
        end do
    else
        gLen=1
    end if

    allocate(gOutH(gLen))

    inquire(IOLENGTH=dbSize) db

    if (myID==rootID)  &
       open(99,file=filename,STATUS='REPLACE',FORM='UNFORMATTED',  &
          ACCESS='direct',RECL=dbSize,IOSTAT=info)

    call MPI_BCast(info,1,MPI_INTEGER,rootID, comm, ierr)

    if (info==0) then 
       if (myID==rootID) then
           ! save data
           disp(1)=1
           do i=1,nNodes-1
             disp(i+1)=disp(i)+locDim(i)
           end do

           dstID=rootID; recInd=disp(dstID+1)
           do i = 1, locDim(dstID+1)
              write(99,rec=recInd) data1(i)
              recInd = recInd+1
           end do

           do i = 2, nNodes
              call MPI_RECV(gOutH,gLen,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,  &
                         MY_DATA, MPI_COMM_WORLD, status, ierr)
              dstID=status(MPI_SOURCE);     recInd = disp(dstID+1)
              do j = 1, locDim(dstID+1)
                 write(99,rec=recInd) gOutH(j)
                 recInd = recInd + 1
              end do
           end do
           close(99)
       else
           call MPI_SEND(data1,N,MPI_DOUBLE_PRECISION,rootID,MY_DATA,  &
                         MPI_COMM_WORLD,ierr)
       end if
  end if
  
  deallocate(gOutH)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PSaveDataSeq_CX(myID,rootID,nNodes,locDim,comm,filename,N,data1,ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer, intent(IN)  :: myID,rootID,nNodes
    integer, intent(IN)  :: locDim(nNodes), comm, N
    double complex, intent(IN) :: data1(N)
    integer, intent(OUT) :: ierr
    integer, parameter :: MY_DATA=1000

    integer :: dbSize, status(MPI_STATUS_SIZE)
    double complex :: db
    double complex, allocatable :: gOutH(:)
    integer :: disp(nNodes), glen, info
    integer :: i,j,dstID, recInd

    if (myID==rootID) then
        glen=locDim(1)
        do i = 2, nNodes
           if (glen>locDim(i)) glen=locDim(i)
        end do
    else
        gLen=1
    end if

    allocate(gOutH(gLen))

    inquire(IOLENGTH=dbSize) db

    if (myID==rootID)  &
       open(99,file=filename,STATUS='REPLACE',FORM='UNFORMATTED',  &
          ACCESS='direct',RECL=dbSize,IOSTAT=info)

    call MPI_BCast(info,1,MPI_INTEGER,rootID, comm, ierr)

    if (info==0) then
       if (myID==rootID) then
           ! save data
           disp(1)=1
           do i=1,nNodes-1
             disp(i+1)=disp(i)+locDim(i)
           end do

           dstID=rootID; recInd=disp(dstID+1)
           do i = 1, locDim(dstID+1)
              write(99,rec=recInd) data1(i)
              recInd = recInd+1
           end do

           do i = 2, nNodes
              call MPI_RECV(gOutH,gLen,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE,  &
                         MY_DATA, MPI_COMM_WORLD, status, ierr)
              dstID=status(MPI_SOURCE);     recInd = disp(dstID+1)
              do j = 1, locDim(dstID+1)
                 write(99,rec=recInd) gOutH(j)
                 recInd = recInd + 1
              end do
           end do
           close(99)
       else
           call MPI_SEND(data1,N,MPI_DOUBLE_COMPLEX,rootID,MY_DATA,  &
                         MPI_COMM_WORLD,ierr)
       end if
  end if

  deallocate(gOutH)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PLoadDataSeq(myID,rootID,nNodes,locDim,comm,filename,N,data1,ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer, intent(IN)  :: myID,rootID,nNodes
    integer, intent(IN)  :: locDim(nNodes), comm, N
    double precision, intent(OUT) :: data1(N)
    integer, intent(OUT) :: ierr
    integer, parameter :: MY_DATA=1000

    integer :: dbSize, status(MPI_STATUS_SIZE)
    double precision :: db
    double precision, allocatable :: gOutH(:)
    integer :: disp(nNodes), glen, info
    integer :: i,j,dstID, recInd

    if (myID==rootID) then
        glen=locDim(1)
        do i = 2, nNodes
           if (glen>locDim(i)) glen=locDim(i)
        end do
    else
        gLen=1
    end if

    allocate(gOutH(gLen))

    inquire(IOLENGTH=dbSize) db

    if (myID==rootID)  &
       open(99,file=filename,STATUS='OLD',FORM='UNFORMATTED',  &
          ACCESS='direct',RECL=dbSize,IOSTAT=info)

    call MPI_BCast(info,1,MPI_INTEGER,rootID, comm, ierr)

    if (info==0) then
       if (myID==rootID) then
           ! load data
           disp(1)=1
           do i=1,nNodes-1
             disp(i+1)=disp(i)+locDim(i)
           end do

           dstID=0; recInd=disp(dstID+1)
           do i = 1, locDim(dstID+1)
              write(99,rec=recInd) gOutH(i)
              recInd = recInd+1
           end do

           do i = 1, nNodes
              recInd=disp(i+1)
              do j = 1, locDim(i+1)
                 read(99,rec=recInd) gOutH(j)
                 recInd = recInd+1
              end do
               
              if (i==rootID) then               
                call MPI_SEND(data1,locDim(i+1),MPI_DOUBLE_PRECISION,i,MY_DATA,  &
                         MPI_COMM_WORLD,ierr)
              else
                data1(1:locDim(i+1))=gOutH(1:locDim(i+1))
              end if
           end do
           close(99)
       else
           call MPI_RECV(data1,N,MPI_DOUBLE_PRECISION,rootID,  &
                         MY_DATA, MPI_COMM_WORLD, status, ierr)
       end if
  end if

  deallocate(gOutH)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PLoadDataSeq_CX(myID,rootID,nNodes,locDim,comm,filename,N,data1,ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer, intent(IN)  :: myID,rootID,nNodes
    integer, intent(IN)  :: locDim(nNodes), comm, N
    double complex, intent(OUT) :: data1(N)
    integer, intent(OUT) :: ierr
    integer, parameter :: MY_DATA=1000

    integer :: dbSize, status(MPI_STATUS_SIZE)
    double complex :: db
    double complex, allocatable :: gOutH(:)
    integer :: disp(nNodes), glen, info
    integer :: i,j,dstID, recInd

    if (myID==rootID) then
        glen=locDim(1)
        do i = 2, nNodes
           if (glen>locDim(i)) glen=locDim(i)
        end do
    else
        gLen=1
    end if

    allocate(gOutH(gLen))

    inquire(IOLENGTH=dbSize) db

    if (myID==rootID)  &
       open(99,file=filename,STATUS='OLD',FORM='UNFORMATTED',  &
          ACCESS='direct',RECL=dbSize,IOSTAT=info)

    call MPI_BCast(info,1,MPI_INTEGER,rootID, comm, ierr)

    if (info==0) then
       if (myID==rootID) then
           ! load data
           disp(1)=1
           do i=1,nNodes-1
             disp(i+1)=disp(i)+locDim(i)
           end do

           dstID=0; recInd=disp(dstID+1)
           do i = 1, locDim(dstID+1)
              write(99,rec=recInd) gOutH(i)
              recInd = recInd+1
           end do

           do i = 1, nNodes
              recInd=disp(i+1)
              do j = 1, locDim(i+1)
                 read(99,rec=recInd) gOutH(j)
                 recInd = recInd+1
              end do

              if (i==rootID) then
                call MPI_SEND(data1,locDim(i+1),MPI_DOUBLE_COMPLEX,i,MY_DATA,  &
                         MPI_COMM_WORLD,ierr)
              else
                data1(1:locDim(i+1))=gOutH(1:locDim(i+1))
              end if
           end do
           close(99)
       else
           call MPI_RECV(data1,N,MPI_DOUBLE_COMPLEX,rootID,  &
                         MY_DATA, MPI_COMM_WORLD, status, ierr)
       end if
  end if

  deallocate(gOutH)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
