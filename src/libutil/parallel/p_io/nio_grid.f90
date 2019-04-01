!ccccccccccccccccccccccccccccccccccccccccccccc
!c       IO subroutines in MPI               c
!ccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine NSaveDataGrid(comm, filename, pos, pSize, N, NG, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN)  :: comm, N, NG, pSize
    double precision, intent(IN) :: data1(N,NG,NG)
    integer, intent(OUT) :: ierr
    
    integer :: dbSize, fh

    call MPI_TYPE_Size(MPI_DOUBLE_PRECISION, dbSize, ierr)

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE,  &
                       MPI_INFO_NULL,fh,ierr)
    call NWriteDataGrid(fh, dbSize, pos, pSize, N, NG, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine NSaveDataGrid_CX(comm, filename, pos, pSize, N, NG, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: comm, N, NG, pSize
    double complex, intent(IN) :: data1(N,NG,NG)
    integer, intent(OUT) :: ierr

    integer :: cxSize, fh

    call MPI_TYPE_Size(MPI_DOUBLE_COMPLEX, cxSize, ierr)

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE,  &
                       MPI_INFO_NULL,fh,ierr)
    call NWriteDataGrid_CX(fh, cxSize, pos, pSize, N, NG, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine NWriteDataGrid(fh, dbSize, pos, pSize, N, nG, data1, ierr)
    implicit none
    include 'mpif.h'
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: fh, dbSize, N, nG, pSize
    double precision, intent(IN) :: data1(N, nG, nG)
    integer, intent(OUT) :: ierr

    integer :: status(MPI_STATUS_SIZE), filetype
    integer(kind=MPI_OFFSET_KIND) :: offset
  
    offset = (pos-1)*dbSize

    call MPI_TYPE_VECTOR(NG*NG,N,pSize,MPI_DOUBLE_PRECISION,fileType,ierr)     
    call MPI_TYPE_COMMIT(filetype,ierr)

    call MPI_FILE_set_View(fh,offset,MPI_DOUBLE_PRECISION,filetype,    &
                 "native",MPI_INFO_NULL, ierr)
    call MPI_FILE_WRITE(fh,data1,N*NG*NG,MPI_DOUBLE_PRECISION,status,ierr)

    call MPI_TYPE_FREE(filetype,ierr)
  
end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine NWriteDataGrid_CX(fh, cxSize, pos, pSize, N, nG, data1, ierr)
    implicit none
    include 'mpif.h'
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: fh, cxSize, N, nG, pSize
    double complex, intent(IN) :: data1(N, nG, nG)
    integer, intent(OUT) :: ierr

    integer :: status(MPI_STATUS_SIZE), filetype
    integer(kind=MPI_OFFSET_KIND) :: offset
  
    offset = (pos-1)*cxSize

    call MPI_TYPE_VECTOR(NG*NG,N,pSize,MPI_DOUBLE_COMPLEX,fileType,ierr)     
    call MPI_TYPE_COMMIT(filetype,ierr)

    call MPI_FILE_set_View(fh,offset,MPI_DOUBLE_COMPLEX,filetype,    &
                 "native",MPI_INFO_NULL, ierr)
    call MPI_FILE_WRITE(fh,data1,N*NG*NG,MPI_DOUBLE_COMPLEX,status,ierr)

    call MPI_TYPE_FREE(filetype,ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!**********************   Load/Read Grid1 Data  ************************

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine NLoadDataGrid(comm, filename, pos, pSize, N, NG, data1, ierr)
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
subroutine NLoadDataGrid_CX(comm, filename, pos, pSize, N, NG, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: comm, N, NG, pSize
    double complex, intent(OUT) :: data1(N,NG,NG)
    integer, intent(OUT) :: ierr

    integer :: cxSize, fh

    call MPI_TYPE_Size(MPI_DOUBLE_COMPLEX, cxSize, ierr)

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    call NReadDataGrid_CX(fh, cxSize, pos, pSize, N, NG, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine NReadDataGrid(fh, dbSize, pos, pSize, N, nG, data1, ierr)
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
subroutine NReadDataGrid_CX(fh, cxSize, pos, pSize, N, nG, data1, ierr)
    implicit none
    include 'mpif.h'
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: fh, cxSize, N, nG, pSize
    double complex, intent(OUT) :: data1(N, nG, nG)
    integer, intent(OUT) :: ierr

    integer :: status(MPI_STATUS_SIZE), filetype
    integer(kind=MPI_OFFSET_KIND) :: offset
  
    offset = (pos-1)*cxSize

    call MPI_TYPE_VECTOR(NG*NG,N,pSize,MPI_DOUBLE_COMPLEX,fileType,ierr)     
    call MPI_TYPE_COMMIT(filetype,ierr)

    call MPI_FILE_set_View(fh,offset,MPI_DOUBLE_COMPLEX,filetype,    &
                 "native",MPI_INFO_NULL, ierr)
    call MPI_FILE_READ(fh,data1,N*NG*NG,MPI_DOUBLE_COMPLEX,status,ierr)

    call MPI_TYPE_FREE(filetype,ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
