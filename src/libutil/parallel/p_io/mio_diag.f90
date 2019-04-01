!ccccccccccccccccccccccccccccccccccccccccccccc
!c       IO subroutines in MPI               c
!ccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MSaveDataDiag(comm, filename, pos, gSize, N, NG, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos, gSize
    integer, intent(IN)  :: comm, N, NG
    double precision, intent(IN) :: data1(N,NG)
    integer, intent(OUT) :: ierr
    
    integer :: dbSize, fh

    call MPI_TYPE_Size(MPI_DOUBLE_PRECISION, dbSize, ierr)

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                 MPI_INFO_NULL,fh,ierr)
    call MWriteDataDiag(fh, dbSize, pos, gSize, N, NG, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MSaveDataDiag_CX(comm, filename, pos, gSize, N, NG, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos, gSize
    integer, intent(IN) :: comm, N, NG
    double complex, intent(IN) :: data1(N,NG)
    integer, intent(OUT) :: ierr

    integer :: cxSize, fh

    call MPI_TYPE_Size(MPI_DOUBLE_COMPLEX, cxSize, ierr)

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE,  &
               MPI_INFO_NULL,fh,ierr)
    call MWriteDataDiag_CX(fh, cxSize, pos, gSize, N, NG, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MWriteDataDiag(fh, dbSize, pos, gSize, N, nG, data1, ierr)
    implicit none
    include 'mpif.h'
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos, gSize
    integer, intent(IN) :: fh, dbSize, N, nG
    double precision, intent(IN) :: data1(N, nG)
    integer, intent(OUT) :: ierr

    integer :: i, status(MPI_STATUS_SIZE)
    integer(kind=MPI_OFFSET_KIND) :: offset
      
    offset = (pos-1)*dbSize

    do i = 1, NG
       call MPI_File_Write_At(fh, offset, data1(1, i), N, &
                   MPI_DOUBLE_PRECISION, status, ierr)
       offset = offset + dbSize*gSize
    end do

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MWriteDataDiag_CX(fh, cxSize, pos, gSize, N, nG, data1, ierr)
    implicit none
    include 'mpif.h'
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos, gSize
    integer, intent(IN) :: fh, cxSize, N, nG
    double complex, intent(IN) :: data1(N, nG)
    integer, intent(OUT) :: ierr

    integer :: i, status(MPI_STATUS_SIZE)
    integer(kind=MPI_OFFSET_KIND) :: offset
      
    offset = (pos-1)*cxSize

    do i = 1, NG
       call MPI_File_Write_At(fh, offset, data1(1, i), N, &
                   MPI_DOUBLE_COMPLEX, status, ierr)
       offset = offset + cxSize*gSize
    end do

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!**********************   Load/Read Grid Data  ************************

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MLoadDataDiag(comm, filename, pos, gSize, N, NG, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos, gSize
    integer, intent(IN)  :: comm, N, NG
    double precision, intent(OUT) :: data1(N,NG)
    integer, intent(OUT) :: ierr
    
    integer :: dbSize, fh

    call MPI_TYPE_Size(MPI_DOUBLE_PRECISION, dbSize, ierr)

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    call MReadDataDiag(fh, dbSize, pos, gSize, N, NG, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MLoadDataDiag_CX(comm, filename, pos, gSize, N, NG, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos, gSize
    integer, intent(IN) :: comm, N, NG
    double complex, intent(OUT) :: data1(N,NG)
    integer, intent(OUT) :: ierr

    integer :: cxSize, fh

    call MPI_TYPE_Size(MPI_DOUBLE_COMPLEX, cxSize, ierr)

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    call MReadDataDiag_CX(fh, cxSize, pos, gSize, N, NG, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MReadDataDiag(fh, dbSize, pos, gSize, N, nG, data1, ierr)
    implicit none
    include 'mpif.h'
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos, gSize
    integer, intent(IN) :: fh, dbSize, N, nG
    double precision, intent(OUT) :: data1(N, nG)
    integer, intent(OUT) :: ierr

    integer :: i, status(MPI_STATUS_SIZE)
    integer(kind=MPI_OFFSET_KIND) :: offset
      
    offset = (pos-1)*dbSize

    do i = 1, NG
       call MPI_File_Read_At(fh, offset, data1(1, i), N, &
                   MPI_DOUBLE_PRECISION, status, ierr)
       offset = offset + dbSize*gSize
    end do

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MReadDataDiag_CX(fh, cxSize, pos, gSize, N, nG, data1, ierr)
    implicit none
    include 'mpif.h'
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos, gSize
    integer, intent(IN) :: fh, cxSize, N, nG
    double complex, intent(OUT) :: data1(N, nG)
    integer, intent(OUT) :: ierr

    integer :: i, status(MPI_STATUS_SIZE)
    integer(kind=MPI_OFFSET_KIND) :: offset
      
    offset = (pos-1)*cxSize

    do i = 1, NG
       call MPI_File_Read_At(fh, offset, data1(1, i), N, &
                   MPI_DOUBLE_COMPLEX, status, ierr)
          offset = offset + cxSize*gSize
    end do

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
