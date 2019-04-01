!ccccccccccccccccccccccccccccccccccccccccccccc
!c       IO subroutines in MPI               c
!ccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PSaveDataSeqAll(myID,rootID,nNodes,locDim,comm,filename,N,data1,ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer, intent(IN)  :: myID,rootID,nNodes
    integer, intent(IN)  :: locDim(nNodes), comm, N
    double precision, intent(IN) :: data1(N)
    integer, intent(OUT) :: ierr
    
    integer :: dbSize
    double precision :: db
    double precision, allocatable :: gOutH(:)
    integer :: disp(nNodes), glen, info 
    integer :: i
   
    if (myID==rootID) then
        glen=sum(locDim(1:nNodes))
    else
        gLen=1
    end if

    allocate(gOutH(gLen))

    inquire(IOLENGTH=dbSize) gOutH

    if (myID==rootID)  &
       open(99,file=filename,STATUS='REPLACE',FORM='UNFORMATTED',  &
          ACCESS='direct',RECL=dbSize,IOSTAT=info)

    call MPI_BCast(info,1,MPI_INTEGER,rootID, comm, ierr)

    if (info==0) then 
       disp(1)=0
       do i=1,nNodes-1
          disp(i+1)=disp(i)+locDim(i)
       end do

       call MPI_GATHERV(data1(1),N,MPI_DOUBLE_PRECISION,gOutH, &
             locDim,disp,MPI_DOUBLE_PRECISION,rootID,comm,ierr) 
       if (myID==rootID)  &
            write(99,rec=1) gOutH
    
       if (myID==rootID)  close(99)      
    end if

    deallocate(gOutH)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PSaveDataSeqAll_CX(myID,rootID,nNodes,locDim,comm,filename,N,data1,ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer, intent(IN)  :: myID,rootID,nNodes
    integer, intent(IN)  :: locDim(nNodes), comm, N
    double complex, intent(IN) :: data1(N)
    integer, intent(OUT) :: ierr
    
    integer :: dbSize
    double complex :: db
    double complex, allocatable :: gOutH(:)
    integer :: disp(nNodes), glen, info 
    integer :: i
   
    if (myID==rootID) then
        glen=sum(locDim(1:nNodes))
    else
        gLen=1
    end if

    allocate(gOutH(gLen))

    inquire(IOLENGTH=dbSize) gOutH

    if (myID==rootID)  &
       open(99,file=filename,STATUS='REPLACE',FORM='UNFORMATTED',  &
          ACCESS='direct',RECL=dbSize,IOSTAT=info)

    call MPI_BCast(info,1,MPI_INTEGER,rootID, comm, ierr)

    if (info==0) then 
       disp(1)=0
       do i=1,nNodes-1
          disp(i+1)=disp(i)+locDim(i)
       end do

       call MPI_GATHERV(data1(1),N,MPI_DOUBLE_COMPLEX,gOutH, &
             locDim,disp,MPI_DOUBLE_COMPLEX,rootID,comm,ierr) 
       if (myID==rootID)  &
            write(99,rec=1) gOutH
    
       if (myID==rootID)  close(99)      
    end if

    deallocate(gOutH)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PLoadDataSeqAll(myID,rootID,nNodes,locDim,comm,filename,N,data1,ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer, intent(IN)  :: myID,rootID,nNodes
    integer, intent(IN)  :: locDim(nNodes), comm, N
    double precision, intent(OUT) :: data1(N)
    integer, intent(OUT) :: ierr
    
    integer :: dbSize
    double precision :: db
    double precision, allocatable :: gOutH(:)
    integer :: disp(nNodes), glen, info 
    integer :: i
   
    if (myID==rootID) then
        glen=sum(locDim(1:nNodes))
    else
        gLen=1
    end if

    allocate(gOutH(gLen))

    inquire(IOLENGTH=dbSize) gOutH

    if (myID==rootID)  &
       open(99,file=filename,STATUS='OLD',FORM='UNFORMATTED',  &
          ACCESS='direct',RECL=dbSize,IOSTAT=info)

    call MPI_BCast(info,1,MPI_INTEGER,rootID, comm, ierr)

    if (info==0) then 
       disp(1)=0
       do i=1,nNodes-1
          disp(i+1)=disp(i)+locDim(i)
       end do

       if (myID==rootID)  &
            read(99,rec=1) gOutH

       call MPI_SCATTERV(gOutH,locDim,disp,MPI_DOUBLE_PRECISION,data1(1), &
                    locDim,MPI_DOUBLE_PRECISION, rootID, comm, ierr  ) 

       if (myID==rootID)  close(99)      

    end if

    deallocate(gOutH)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PLoadDataSeqAll_CX(myID,rootID,nNodes,locDim,comm,filename,N,data1,ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer, intent(IN)  :: myID,rootID,nNodes
    integer, intent(IN)  :: locDim(nNodes), comm, N
    double complex, intent(OUT) :: data1(N)
    integer, intent(OUT) :: ierr
    
    integer :: dbSize
    double complex :: db
    double complex, allocatable :: gOutH(:)
    integer :: disp(nNodes), glen, info 
    integer :: i
   
    if (myID==rootID) then
        glen=sum(locDim(1:nNodes))
    else
        gLen=1
    end if

    allocate(gOutH(gLen))

    inquire(IOLENGTH=dbSize) gOutH

    if (myID==rootID)  &
       open(99,file=filename,STATUS='OLD',FORM='UNFORMATTED',  &
          ACCESS='direct',RECL=dbSize,IOSTAT=info)

    call MPI_BCast(info,1,MPI_INTEGER,rootID, comm, ierr)

    if (info==0) then 
       disp(1)=0
       do i=1,nNodes-1
          disp(i+1)=disp(i)+locDim(i)
       end do

       if (myID==rootID)  &
            read(99,rec=1) gOutH

       call MPI_SCATTERV(gOutH,locDim,disp,MPI_DOUBLE_COMPLEX,data1(1), &
                    locDim,MPI_DOUBLE_COMPLEX, rootID,comm,ierr  ) 

       if (myID==rootID)  close(99)      

    end if

    deallocate(gOutH)

end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
