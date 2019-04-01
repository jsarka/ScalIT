!ccccccccccccccccccccccccccccccccccccccccccccc
!c       IO subroutines in MPI               c
!ccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MSaveDataSeq(comm, filename, pos, pLen, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN)  :: comm, pLen
    double precision, intent(IN) :: data1(pLen)
    integer, intent(OUT) :: ierr
    
    integer :: fh

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                 MPI_INFO_NULL,fh,ierr)
    call MWriteDataSeq(fh, pos, pLen, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MWriteDataSeq(fh, pos, pLen, data1, ierr)
    implicit none
    include 'mpif.h'
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: fh, pLen
    double precision, intent(IN) :: data1(pLen)
    integer, intent(OUT) :: ierr

    integer :: dbSize, status(MPI_STATUS_SIZE)
    integer(kind=MPI_OFFSET_KIND) :: offset

    call MPI_TYPE_Size(MPI_DOUBLE_PRECISION, dbSize, ierr)
 
    offset = (pos-1)*dbSize
    
    call MPI_File_Write_At(fh, offset, data1(1), pLen, &
                   MPI_DOUBLE_PRECISION, status, ierr)   

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MSaveDataSeq_CX(comm, filename, pos, pLen, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: comm, plen
    double complex, intent(IN) :: data1(plen)
    integer, intent(OUT) :: ierr

    integer :: fh

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE,  &
               MPI_INFO_NULL,fh,ierr)
    call MWriteDataSeq_CX(fh, pos, pLen, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MWriteDataSeq_CX(fh, pos, pLen, data1, ierr)
    implicit none
    include 'mpif.h'
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: fh, pLen
    double complex, intent(IN) :: data1(plen)
    integer, intent(OUT) :: ierr

    integer :: cxSize, status(MPI_STATUS_SIZE)
    integer(kind=MPI_OFFSET_KIND) :: offset
      
    call MPI_TYPE_Size(MPI_DOUBLE_COMPLEX, cxSize, ierr)

    offset = (pos-1)*cxSize

    call MPI_File_Write_At(fh, offset, data1(1), pLen, &
                   MPI_DOUBLE_COMPLEX, status, ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!**********************   Load/Read Grid Data  ************************

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MLoadDataSeq(comm, filename, pos, pLen, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN)  :: comm,pLen
    double precision, intent(OUT) :: data1(pLen)
    integer, intent(OUT) :: ierr
    
    integer :: fh

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    call MReadDataSeq(fh, pos, pLen, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MReadDataSeq(fh, pos, pLen, data1, ierr)
    implicit none
    include 'mpif.h'
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: fh, pLen
    double precision, intent(OUT) :: data1(pLen)
    integer, intent(OUT) :: ierr

    integer ::  dbSize, status(MPI_STATUS_SIZE)
    integer(kind=MPI_OFFSET_KIND) :: offset      

    call MPI_TYPE_Size(MPI_DOUBLE_PRECISION, dbSize, ierr)

    offset = (pos-1)*dbSize

    call MPI_File_Read_At(fh, offset, data1(1), pLen, &
                   MPI_DOUBLE_PRECISION, status, ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MLoadDataSeq_CX(comm, filename, pos, pLen, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: comm, pLen
    double complex, intent(OUT) :: data1(pLen)
    integer, intent(OUT) :: ierr

    integer :: fh

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    call MReadDataSeq_CX(fh, pos, pLen, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MReadDataSeq_CX(fh, pos, pLen, data1, ierr)
    implicit none
    include 'mpif.h'
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: fh,  pLen
    double complex, intent(OUT) :: data1(pLen)
    integer, intent(OUT) :: ierr

    integer ::  cxSize, status(MPI_STATUS_SIZE)
    integer(kind=MPI_OFFSET_KIND) :: offset

    call MPI_TYPE_Size(MPI_DOUBLE_COMPLEX, cxSize, ierr)
 
    offset = (pos-1)*cxSize

    call MPI_File_Read_At(fh, offset, data1(1), plen, &
                   MPI_DOUBLE_COMPLEX, status, ierr)   

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
