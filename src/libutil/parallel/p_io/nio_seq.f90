!ccccccccccccccccccccccccccccccccccccccccccccc
!c       IO subroutines in MPI               c
!ccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine NSaveDataSeq(comm, filename, pos, pSize, NG, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN)  :: comm, NG, pSize
    double precision, intent(IN) :: data1(NG)
    integer, intent(OUT) :: ierr
    
    integer :: dbSize, fh

    call MPI_TYPE_Size(MPI_DOUBLE_PRECISION, dbSize, ierr)

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE,  &
                       MPI_INFO_NULL,fh,ierr)
    call NWriteDataSeq(fh, dbSize, pos, pSize, NG, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine NSaveDataSeq_CX(comm, filename, pos, pSize, NG, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: comm, NG, pSize
    double complex, intent(IN) :: data1(NG)
    integer, intent(OUT) :: ierr

    integer :: cxSize, fh

    call MPI_TYPE_Size(MPI_DOUBLE_COMPLEX, cxSize, ierr)

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                       MPI_INFO_NULL,fh,ierr)
    call NWriteDataSeq_CX(fh, cxSize, pos, pSize, NG, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine NWriteDataSeq(fh, dbSize, pos, pSize, nG, data1, ierr)
    implicit none
    include 'mpif.h'
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: fh, dbSize, nG, pSize
    double precision, intent(IN) :: data1( nG)
    integer, intent(OUT) :: ierr

    integer :: status(MPI_STATUS_SIZE), filetype
    integer(kind=MPI_OFFSET_KIND) :: offset
  
    offset = (pos-1)*dbSize

    call MPI_TYPE_VECTOR(NG,1,pSize,MPI_DOUBLE_PRECISION,fileType,ierr)     
    call MPI_TYPE_COMMIT(filetype,ierr)

    call MPI_FILE_set_View(fh,offset,MPI_DOUBLE_PRECISION,filetype,    &
                 "native",MPI_INFO_NULL, ierr)
    call MPI_FILE_WRITE(fh,data1,NG,MPI_DOUBLE_PRECISION,status,ierr)

    call MPI_TYPE_FREE(filetype,ierr)
  
end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine NWriteDataSeq_CX(fh, cxSize, pos, pSize, nG, data1, ierr)
    implicit none
    include 'mpif.h'
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: fh, cxSize, nG, pSize
    double complex, intent(IN) :: data1(nG)
    integer, intent(OUT) :: ierr

    integer :: status(MPI_STATUS_SIZE), filetype
    integer(kind=MPI_OFFSET_KIND) :: offset
  
    offset = (pos-1)*cxSize

    call MPI_TYPE_VECTOR(NG,1,pSize,MPI_DOUBLE_COMPLEX,fileType,ierr)     
    call MPI_TYPE_COMMIT(filetype,ierr)

    call MPI_FILE_set_View(fh,offset,MPI_DOUBLE_COMPLEX,filetype,    &
                 "native",MPI_INFO_NULL, ierr)
    call MPI_FILE_WRITE(fh,data1,NG,MPI_DOUBLE_COMPLEX,status,ierr)

    call MPI_TYPE_FREE(filetype,ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!**********************   Load/Read diag1 Data  ************************

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine NLoadDataSeq(comm, filename, pos, pSize, NG, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN)  :: comm, NG, pSize
    double precision, intent(OUT) :: data1(NG)
    integer, intent(OUT) :: ierr
    
    integer :: dbSize, fh

    call MPI_TYPE_Size(MPI_DOUBLE_PRECISION, dbSize, ierr)

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    call NReadDataSeq(fh, dbSize, pos, pSize, NG, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine NLoadDataSeq_CX(comm, filename, pos, pSize, NG, data1, ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: comm, NG, pSize
    double complex, intent(OUT) :: data1(NG)
    integer, intent(OUT) :: ierr

    integer :: cxSize, fh

    call MPI_TYPE_Size(MPI_DOUBLE_COMPLEX, cxSize, ierr)

    call MPI_FILE_OPEN(comm,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    call NReadDataSeq_CX(fh, cxSize, pos, pSize, NG, data1, ierr)
    call MPI_FILE_CLOSE(fh, ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine NReadDataSeq(fh, dbSize, pos, pSize, nG, data1, ierr)
    implicit none
    include 'mpif.h'
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: fh, dbSize, nG, pSize
    double precision, intent(OUT) :: data1(nG)
    integer, intent(OUT) :: ierr

    integer :: status(MPI_STATUS_SIZE), filetype
    integer(kind=MPI_OFFSET_KIND) :: offset
  
    offset = (pos-1)*dbSize

    call MPI_TYPE_VECTOR(NG,1,pSize,MPI_DOUBLE_PRECISION,fileType,ierr)     
    call MPI_TYPE_COMMIT(filetype,ierr)

    call MPI_FILE_set_View(fh,offset,MPI_DOUBLE_PRECISION,filetype,    &
                 "native",MPI_INFO_NULL, ierr)
    call MPI_FILE_READ(fh,data1,NG,MPI_DOUBLE_PRECISION,status,ierr)

    call MPI_TYPE_FREE(filetype,ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine NReadDataSeq_CX(fh, cxSize, pos, pSize, nG, data1, ierr)
    implicit none
    include 'mpif.h'
    integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos
    integer, intent(IN) :: fh, cxSize, nG, pSize
    double complex, intent(OUT) :: data1(nG)
    integer, intent(OUT) :: ierr

    integer :: status(MPI_STATUS_SIZE), filetype
    integer(kind=MPI_OFFSET_KIND) :: offset
  
    offset = (pos-1)*cxSize

    call MPI_TYPE_VECTOR(NG,1,pSize,MPI_DOUBLE_COMPLEX,fileType,ierr)     
    call MPI_TYPE_COMMIT(filetype,ierr)

    call MPI_FILE_set_View(fh,offset,MPI_DOUBLE_COMPLEX,filetype,    &
                 "native",MPI_INFO_NULL, ierr)
    call MPI_FILE_READ(fh,data1,NG,MPI_DOUBLE_COMPLEX,status,ierr)

    call MPI_TYPE_FREE(filetype,ierr)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
