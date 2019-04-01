!ccccccccccccccccccccccccccccccccccccccccccccc
!c       IO subroutines in MPI               c
!ccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PMSaveDataGrid(myID, rootID, nNodes, comm, filename, pos, gSize, N, NG, data1, ierr)
    implicit none
    include 'mpif.h'
    integer, intent(IN) :: myID,rootID,nNodes
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos, gSize
    integer, intent(IN)  :: comm, N, NG
    double precision, intent(IN) :: data1(N,NG,NG)
    integer, intent(OUT) :: ierr
    
    integer :: dbSize, fh

    call MPI_TYPE_Size(MPI_DOUBLE_PRECISION, dbSize, ierr)

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                       MPI_INFO_NULL,fh,ierr)
    call PMWriteDataGrid(myID, rootID, nNodes, fh, dbSize, pos, gSize, N, NG, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PMWriteDataGrid(myID, rootID, nNodes, fh, dbSize, pos, gSize, N, nG, data1, ierr)
    implicit none
    include 'mpif.h'
    integer, intent(IN) :: myID,rootID,nNodes
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos, gSize
    integer, intent(IN) :: fh, dbSize, N, nG
    double precision, intent(IN) :: data1(N, nG, nG)
    integer, intent(OUT) :: ierr

    integer :: i, j, k, status(MPI_STATUS_SIZE)
    integer(kind=MPI_OFFSET_KIND) :: offset

    if (myid==rootID) then
       print *, 'PMWriteDataGrid is running!'
       print *
    end if

    offset = (pos-1)*dbSize

    do k = 0, nNodes-1

       if (myID==k) then
          print *, k, 'th tread is writing'

          do i = 1, NG
             do j = 1, NG
                call MPI_File_Write_At(fh, offset, data1(1, j, i), N, &
                         MPI_DOUBLE_PRECISION, status, ierr)
                offset = offset + dbSize*gSize
             end do
          end do
       end if
       call MPI_Barrier(MPI_COMM_WORLD, ierr)

    end do


end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!**********************   Load/Read Grid Data  ************************

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PMLoadDataGrid(comm, filename, pos, gSize, N, NG, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos, gSize
    integer, intent(IN)  :: comm, N, NG
    double precision, intent(OUT) :: data1(N,NG,NG)
    integer, intent(OUT) :: ierr
    
    integer :: dbSize, fh
    
    call MPI_TYPE_Size(MPI_DOUBLE_PRECISION, dbSize, ierr)

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    call PMReadDataGrid(fh, dbSize, pos, gSize, N, NG, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PMReadDataGrid(fh, dbSize, pos, gSize, N, nG, data1, ierr)
    implicit none
    include 'mpif.h'
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos, gSize
    integer, intent(IN) :: fh, dbSize, N, nG
    double precision, intent(OUT) :: data1(N, nG, nG)
    integer, intent(OUT) :: ierr

    integer :: i, j, status(MPI_STATUS_SIZE)
    integer(kind=MPI_OFFSET_KIND) :: offset
      
    offset = (pos-1)*dbSize

    do i = 1, NG
       do j = 1, NG
          call MPI_File_Read_At(fh, offset, data1(1, j, i), N, &
                   MPI_DOUBLE_PRECISION, status, ierr)
          offset = offset + dbSize*gSize
       end do
    end do

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
