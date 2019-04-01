!ccccccccccccccccccccccccccccccccccccccccccccc
!c       IO subroutines in MPI               c
!ccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PNSaveDataGrid(myID, rootID, nNodes, comm, filename, pos, pSize, N, NG, data1, ierr)
    implicit none
    include 'mpif.h'
    integer, intent(IN) :: myID,rootID,nNodes
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN)  :: comm, N, NG, pSize
    double precision, intent(IN) :: data1(N,NG,NG)
    integer, intent(OUT) :: ierr
    
    integer :: dbSize, fh

    call MPI_TYPE_Size(MPI_DOUBLE_PRECISION, dbSize, ierr)

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE,  &
                       MPI_INFO_NULL,fh,ierr)
    call PNWriteDataGrid(myID, rootID, nNodes, fh, dbSize, pos, pSize, N, NG, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PNWriteDataGrid(myID, rootID, nNodes, fh, dbSize, pos, pSize, N, nG, data1, ierr)
    implicit none
    include 'mpif.h'
    integer, intent(IN) :: myID,rootID,nNodes
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: fh, dbSize, N, nG, pSize
    double precision, intent(IN) :: data1(N, nG, nG)
    integer, intent(OUT) :: ierr

    integer :: i, j, k, status(MPI_STATUS_SIZE), filetype
    integer(kind=MPI_OFFSET_KIND) :: offset

    if (myid==rootID) then
       print *, 'PNWriteDataGrid is running!'
       print *
    end if

    offset = (pos-1)*dbSize

    call MPI_TYPE_VECTOR(NG*NG,N,pSize,MPI_DOUBLE_PRECISION,fileType,ierr)     
    call MPI_TYPE_COMMIT(filetype,ierr)

    call MPI_FILE_set_View(fh,offset,MPI_DOUBLE_PRECISION,filetype,    &
                 "native",MPI_INFO_NULL, ierr)

    do k = 0, nNodes-1
       if (myID==k) then
          print *, myid, 'th tread is writing'
          call MPI_FILE_WRITE(fh,data1,N*NG*NG,MPI_DOUBLE_PRECISION,status,ierr)
       end if
       call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end do

    call MPI_TYPE_FREE(filetype,ierr)
  
end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!**********************   Load/Read Grid1 Data  ************************

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PNLoadDataGrid(comm, filename, pos, pSize, N, NG, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN)  :: comm, N, NG, pSize
    double precision, intent(OUT) :: data1(N,NG,NG)
    integer, intent(OUT) :: ierr
    
    integer :: dbSize, fh

    call MPI_TYPE_Size(MPI_DOUBLE_PRECISION, dbSize, ierr)

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    call NReadDataGrid(fh, dbSize, pos, pSize, N, NG, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PNReadDataGrid(fh, dbSize, pos, pSize, N, nG, data1, ierr)
    implicit none
    include 'mpif.h'
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: fh, dbSize, N, nG, pSize
    double precision, intent(OUT) :: data1(N, nG, nG)
    integer, intent(OUT) :: ierr

    integer :: status(MPI_STATUS_SIZE), filetype
    integer(kind=MPI_OFFSET_KIND) :: offset
  
    offset = (pos-1)*dbSize

    call MPI_TYPE_VECTOR(NG*NG,N,pSize,MPI_DOUBLE_PRECISION,fileType,ierr)     
    call MPI_TYPE_COMMIT(filetype,ierr)

    call MPI_FILE_set_View(fh,offset,MPI_DOUBLE_PRECISION,filetype,    &
                 "native",MPI_INFO_NULL, ierr)
    call MPI_FILE_READ(fh,data1,N*NG*NG,MPI_DOUBLE_PRECISION,status,ierr)

    call MPI_TYPE_FREE(filetype,ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
