!ccccccccccccccccccccccccccccccccccccccccccccc
!c       IO subroutines in MPI               c
!ccccccccccccccccccccccccccccccccccccccccccccc

!******************   Save/Write data to MPI File **********************

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MSaveData(comm, filename, pos, N, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: comm, N
    double precision, intent(IN) :: data1(N)
    integer, intent(OUT) :: ierr

    integer :: fh, dbSize

    call MPI_TYPE_Size(MPI_DOUBLE_PRECISION, dbSize, ierr)

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                       MPI_INFO_NULL,fh,ierr)
    call MWriteData(fh, dbSize, pos, N, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MSaveData_CX(comm, filename,pos, N, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: comm, N
    double complex, intent(IN) :: data1(N)
    integer, intent(OUT) :: ierr

    integer :: fh, cxSize

    call MPI_TYPE_Size(MPI_DOUBLE_COMPLEX, cxSize, ierr)

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE,    &
                      MPI_INFO_NULL,fh,ierr)
    call MWriteData_CX(fh, cxSize, pos,  N, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MWriteData(fh, dbSize, pos,  N, data1, ierr)
    implicit none
    include 'mpif.h'    
    integer, intent(IN) :: N, fh, dbSize
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    double precision, intent(IN) :: data1(N)
    integer, intent(OUT) :: ierr

    integer :: status(MPI_STATUS_SIZE)
    integer(kind=MPI_OFFSET_KIND) :: offset

    offset = (pos-1)*dbSize

    call MPI_FILE_WRITE_AT(fh,offset,data1,N,MPI_DOUBLE_PRECISION,status,ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MWriteData_CX(fh, cxSize, pos, N, data1, ierr)
    implicit none
    include 'mpif.h'    
    integer, intent(IN) :: N, fh, cxSize
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    double complex, intent(IN) :: data1(N)
    integer, intent(OUT) :: ierr

    integer :: status(MPI_STATUS_SIZE)
    integer(kind=MPI_OFFSET_KIND) :: offset

    offset = (pos-1)*cxSize

    call MPI_FILE_WRITE_AT(fh,offset,data1,N,MPI_DOUBLE_COMPLEX,status,ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



!************************  Load/Read data from MPI file  ***************
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MLoadData(comm, filename, pos, N, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: comm, N
    double precision, intent(OUT) :: data1(N)
    integer, intent(OUT) :: ierr

    integer :: fh, dbSize

    call MPI_TYPE_Size(MPI_DOUBLE_PRECISION, dbSize, ierr)

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    call MReadData(fh, dbSize, pos, N, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MLoadData_CX(comm, filename,pos, N, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: comm, N
    double complex, intent(OUT) :: data1(N)
    integer, intent(OUT) :: ierr

    integer :: fh, cxSize

    call MPI_TYPE_Size(MPI_DOUBLE_COMPLEX, cxSize, ierr)

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    call MReadData_CX(fh, cxSize, pos, N, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MLoadData1(filename, pos, N, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: N
    double precision, intent(OUT) :: data1(N)
    integer, intent(OUT) :: ierr

    integer :: fh, dbSize

    call MPI_TYPE_Size(MPI_DOUBLE_PRECISION, dbSize, ierr)

    call MPI_FILE_OPEN(MPI_COMM_WORLD,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    call MReadData(fh, dbSize, pos,  N, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MLoadData1_CX(filename,pos, N, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: N
    double complex, intent(OUT) :: data1(N)
    integer, intent(OUT) :: ierr

    integer :: fh, cxSize

    call MPI_TYPE_Size(MPI_DOUBLE_COMPLEX, cxSize, ierr)

    call MPI_FILE_OPEN(MPI_COMM_WORLD,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    call MReadData_CX(fh, cxSize,  pos,  N, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MReadData(fh, dbSize, pos, N, data1, ierr)
    implicit none
    include 'mpif.h'    
    integer, intent(IN) :: N, fh, dbSize
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    double precision, intent(OUT) :: data1(N)
    integer, intent(OUT) :: ierr

    integer :: status(MPI_STATUS_SIZE)
    integer(kind=MPI_OFFSET_KIND) :: offset

    offset = (pos-1)*dbSize

    call MPI_FILE_READ_AT(fh,offset,data1,N,MPI_DOUBLE_PRECISION,status,ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MReadData_CX(fh, cxSize, pos, N, data1, ierr)
    implicit none
    include 'mpif.h'    
    integer, intent(IN) :: N, fh, cxSize
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    double complex, intent(OUT) :: data1(N)
    integer, intent(OUT) :: ierr

    integer :: status(MPI_STATUS_SIZE)
    integer(kind=MPI_OFFSET_KIND) :: offset

    offset = (pos-1)*cxSize

    call MPI_FILE_READ_AT(fh,offset,data1,N,MPI_DOUBLE_COMPLEX,status,ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

