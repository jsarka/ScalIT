!ccccccccccccccccccccccccccccccccccccccccccccc
!c       IO subroutines in MPI               c
!ccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PSaveDataGrid(myID,rootID,nNodes,locDim,comm,filename,N,NG,data1,ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer, intent(IN) :: myID,rootID,nNodes
    integer, intent(IN) :: locDim(nNodes), comm, N, NG
    double precision, intent(IN) :: data1(N,NG,NG)
    integer, intent(OUT) :: ierr
    
    integer :: dbSize
    double precision :: db
    double precision, allocatable :: gOutH(:)
    integer :: disp(nNodes)
    integer :: glen, recLen, info, recInd, locLen      
    integer :: i, j

    glen=sum(locDim(1:nNodes))
    locLen=locDim(myId+1)
    if (rootID==myID) then
      allocate(gOutH(glen))
    else
      allocate(gOutH(1))
    end if

   inquire(IOLENGTH=recLen) gOutH
   if (myID==rootID)  &
     open(99,file=filename,STATUS='Replace',FORM='UNFORMATTED',  &
          ACCESS='direct',RECL=recLen,IOSTAT=info)

   call MPI_BCAST(info,1,MPI_INTEGER,rootID,comm,ierr)

   if (info==0) then 
      disp(1)=0
      do i=1,nNodes-1
          disp(i+1)=disp(i)+locDim(i)
      end do

      recInd=1       
      do i = 1, NG
         do j = 1, NG
            call MPI_GATHERV(data1(1,j,i),N,MPI_DOUBLE_PRECISION,gOutH, &
                  locDim,disp,MPI_DOUBLE_PRECISION,rootID,comm,ierr) 
            if (myID==rootID) write(99,rec=recInd) gOutH
            recInd=recInd+1
         end do
      end do
    
      if (myID==rootID)  close(99)      
   end if

   if (allocated(gOutH))  deallocate(gOutH)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PSaveDataGrid_CX(myID,rootID,nNodes,locDim,comm,filename,N,NG,data1,ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer, intent(IN) :: myID,rootID,nNodes
    integer, intent(IN) :: locDim(nNodes), comm, N, NG
    double complex, intent(IN) :: data1(N,NG,NG)
    integer, intent(OUT) :: ierr
    
    integer :: dbSize
    double complex :: db
    double complex, allocatable :: gOutH(:)
    integer :: disp(nNodes)
    integer :: glen, recLen, info, recInd, locLen       
    integer :: i, j

    glen=sum(locDim(1:nNodes))
    locLen=locDim(myId+1)
    if (rootID==myID) then
      allocate(gOutH(glen))
    else
      allocate(gOutH(1))
    end if

   inquire(IOLENGTH=recLen) gOutH
   if (myID==rootID)  &
     open(99,file=filename,STATUS='Replace',FORM='UNFORMATTED',  &
          ACCESS='direct',RECL=recLen,IOSTAT=info)

   call MPI_BCAST(info,1,MPI_INTEGER,rootID,comm, ierr)

   if (info==0) then 
      disp(1)=0
      do i=1,nNodes-1
          disp(i+1)=disp(i)+locDim(i)
      end do

      recInd=1       
      do i = 1, NG
         do j = 1, NG
            call MPI_GATHERV(data1(1,j,i),N,MPI_DOUBLE_COMPLEX,gOutH, &
                  locDim,disp,MPI_DOUBLE_COMPLEX,rootID,comm,ierr) 
            if (myID==rootID)  &
                write(99,rec=recInd) gOutH
            recInd=recInd+1
         end do
      end do
    
      if (myID==rootID)  close(99)      
   end if

   if (allocated(gOutH))  deallocate(gOutH)

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!**********************   Load/Read Grid Data  ************************

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PLoadDataGrid(myID,rootID,nNodes,locDim,comm,filename,N,NG,data1,ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer, intent(IN) :: myID,rootID,nNodes
    integer, intent(IN) :: locDim(nNodes), comm, N, NG
    double precision, intent(OUT) :: data1(N,NG,NG)
    integer, intent(OUT) :: ierr
    
    integer :: dbSize
    double precision :: db
    double precision, allocatable :: gOutH(:)
    integer :: disp(nNodes)
    integer :: glen, recLen, info, recInd, locLen       
    integer :: i, j

    glen=sum(locDim(1:nNodes))
    locLen=locDim(myId+1)
    if (rootID==myID) then
      allocate(gOutH(glen))
    else
      allocate(gOutH(1))
    end if

   inquire(IOLENGTH=recLen) gOutH
   if (myID==rootID)  &
     open(99,file=filename,STATUS='OLD',FORM='UNFORMATTED',  &
          ACCESS='direct',RECL=recLen,IOSTAT=info)

   call MPI_BCAST(info,1,MPI_INTEGER,rootID,comm,ierr)

   if (info==0) then 
      disp(1)=0
      do i=1,nNodes-1
          disp(i+1)=disp(i)+locDim(i)
      end do

      recInd=1       
      do i = 1, NG
         do j = 1, NG
            if (myID==rootID)  &
                read(99,rec=recInd) gOutH
            call MPI_SCATTERV(gOutH,locDim,disp,MPI_DOUBLE_PRECISION,data1(1,j,i), &
                    locLen,MPI_DOUBLE_PRECISION,rootID,comm,ierr) 
            recInd=recInd+1
         end do
      end do
    
      if (myID==rootID)  close(99)      
   end if

   if (allocated(gOutH))  deallocate(gOutH)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine PLoadDataGrid_CX(myID,rootID,nNodes,locDim,comm,filename,N,NG,data1,ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer, intent(IN) :: myID,rootID,nNodes
    integer, intent(IN) :: locDim(nNodes), comm, N, NG
    double complex, intent(OUT) :: data1(N,NG,NG)
    integer, intent(OUT) :: ierr
    
    integer :: dbSize
    double complex :: db
    double complex, allocatable :: gOutH(:)
    integer :: disp(nNodes)
    integer :: glen, recLen, info, recInd, locLen       
    integer :: i, j

    glen=sum(locDim(1:nNodes))
    locLen=locDim(myId+1)
    if (rootID==myID) then
      allocate(gOutH(glen))
    else
      allocate(gOutH(1))
    end if

   inquire(IOLENGTH=recLen) gOutH
   if (myID==rootID)  &
     open(99,file=filename,STATUS='OLD',FORM='UNFORMATTED',  &
          ACCESS='direct',RECL=recLen,IOSTAT=info)

   call MPI_BCAST(info,1,MPI_INTEGER,rootID,comm,ierr)

   if (info==0) then 
      disp(1)=0
      do i=1,nNodes-1
          disp(i+1)=disp(i)+locDim(i)
      end do

      recInd=1       
      do i = 1, NG
         do j = 1, NG
            if (myID==rootID)  &
                read(99,rec=recInd) gOutH
            call MPI_SCATTERV(gOutH,locDim,disp,MPI_DOUBLE_COMPLEX,data1(1,j,i), &
                    locLen,MPI_DOUBLE_COMPLEX,rootID,comm,ierr) 
            recInd=recInd+1
         end do
      end do
    
      if (myID==rootID)  close(99)      
   end if

   if (allocated(gOutH))  deallocate(gOutH)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
