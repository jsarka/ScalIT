!ccccccccccccccccccccccccccccccccccccccccccccc
!c       IO subroutines in MPI               c
!ccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine P2SaveDataGrid(myID,rootID,nNodes,locDim,comm,filename,N,NG,data1,ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer, intent(IN) :: myID,rootID,nNodes
    integer, intent(IN) :: locDim(nNodes), comm, N, NG
    double precision, intent(IN) :: data1(N,NG,NG)
    integer, intent(OUT) :: ierr
    
    integer :: dbSize
    double precision :: db
    double precision, allocatable :: gOutH(:), data0(:)
    integer :: disp(nNodes)
    integer :: glen, recLen, info, recInd, locLen      
    integer :: locStart(nNodes)

    integer :: i, j, k
    character(len=8) :: fs, is

    character(len=20) :: stmp
    character(len=40) :: filename2

    if (myid==rootID) then
       print *, 'P2SaveDataGrid is running!'
       print *
    end if

    glen=sum(locDim(1:nNodes))
    locLen=locDim(myId+1)

    locStart(1)=1;
    do i = 1, nNodes-1
       locStart(i+1)=locStart(i)+locDim(i)
    end do

   if (myID==rootID) then 
     open (unit=999,file='input/param_h0gm.txt')

     write(999,*) "rNum, jkNum, nproc"
     write(999,*)  glen, NG, nNodes
     write(999,*)

     write(999,*) 'locDim'
     write(999,*) locDim
     !do j = 1, nNodes
     !   write(999,*) locDim(j)
     !enddo
     !write(999,*)

      close(999)
   end if 

   !fs = "(I8)"
   !write (is,fs) myid
   !is=adjustl(is)
   !filename2 = filename//'_'//trim(is)

   write (stmp, "(A1,I0.5,A4)") "_", myid+1, ".dat"
   filename2 = trim(filename)//trim(stmp)
   filename2 = trim(filename2)

   do k = 0, nNodes-1

     if (myID==k) then
       print *, k, 'th tread is writing'
       !allocate(data0(NG,NG))
       inquire(IOLENGTH=recLen) data1(1,:,:)

       open(k+99,file=filename,STATUS='UNKNOWN',FORM='UNFORMATTED',  &
            ACCESS='direct',RECL=recLen,IOSTAT=info)

         if (info==0) then 
            recInd=locStart(k+1)
            do j = 1, locDim(myID+1)
               write(k+99,rec=recInd) data1(j,1:NG,1:NG)
               recInd=recInd+1
            end do
            !write(99,rec=1) data1
            close(k+99)      
         end if

         !deallocate(data0)

     end if ! myID==k
     call MPI_Barrier(MPI_COMM_WORLD, ierr)

   end do

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine P2SaveDataGrid_CX(myID,rootID,nNodes,locDim,comm,filename,N,NG,data1,ierr)
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
subroutine P2LoadDataGrid(myID,rootID,nNodes,locDim,comm,filename,N,NG,data1,ierr)
    implicit none
    include 'mpif.h'
    character(len=*), intent(IN) :: filename
    integer, intent(IN) :: myID,rootID,nNodes
    integer, intent(IN) :: locDim(nNodes), comm, N, NG
    double precision, intent(OUT) :: data1(N,NG,NG)
    integer, intent(OUT) :: ierr
    
    integer :: dbSize
    double precision :: db
    double precision, allocatable :: gOutH(:), data0(:)
    integer :: disp(nNodes)
    integer :: glen, recLen, info, recInd, locLen       
    integer :: i, j, k

    glen=sum(locDim(1:nNodes))
    locLen=locDim(myId+1)

    disp(1)=0
    do i=1,nNodes-1
      disp(i+1)=disp(i)+locDim(i)
    end do

   do k = 0, nNodes-1

     if (myID==k) then
       print *, k, 'th tread is reading'
       inquire(IOLENGTH=recLen) data1(1,:,:)

       open(k+99,file=filename,STATUS='OLD',FORM='UNFORMATTED',  &
            ACCESS='direct',RECL=recLen,IOSTAT=info)

         if (info==0) then
            recInd = disp(k+1)+1
            do j = 1, locDim(myID+1)
               read(k+99,rec=recInd) data1(j,1:NG,1:NG)
               recInd=recInd+1
            end do
            !read(99,rec=1) data1
            close(99)
         end if

         !deallocate(data0)

     end if ! myID==k
     call MPI_Barrier(MPI_COMM_WORLD, ierr)

   end do

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine P2LoadDataGrid_CX(myID,rootID,nNodes,locDim,comm,filename,N,NG,data1,ierr)
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
