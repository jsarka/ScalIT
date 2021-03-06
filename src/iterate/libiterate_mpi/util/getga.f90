!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                 Initialize the global array: x0(i)=i                c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine InitGA(level, sF, sN, nNodes, id, Nin, Nout, locData )
    implicit none
    include 'mpif.h'
    integer, intent(IN) :: level, sF, sN(sF), nNodes, id, Nin, Nout
    double precision, intent(OUT) :: locData(Nin, nout)

    integer(kind=MPI_OFFSET_KIND) :: sLen,sPos(nNodes),ePos(nNodes)
    integer :: bNum(nNodes), locDim(nNodes)
    integer :: i, j, sff, soff

    call getPosition(level, sF, sN, nNodes, sLen, sPos, ePos, bNum, locDim)    

    do i = 1, Nout
       sff = sPos(id+1) + (i-1)*sLen
       do j = 1, Nin
          locData(j,i) = sff + (j-1)
       end do
    end do
 
end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                  Get the content of global array                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getGA(level, sF, sN, nNodes, id, N1, locData, N, gaData )
    implicit none
    include 'mpif.h'
    integer, intent(IN) :: level, sF, sN(sF), nNodes, id,N1, N
    double precision, intent(IN)  :: locData(N1)
    double precision, intent(OUT) :: gaData(N)

    integer, parameter :: DAT_TAG = 200
    integer(kind=MPI_OFFSET_KIND) :: sLen,sPos(nNodes),ePos(nNodes)
    integer :: bNum(nNodes), locDim(nNodes)
    integer :: ierr, status(MPI_STATUS_SIZE)
    integer :: i, j, pos,pos1, len

    call getPosition(level, sF, sN, nNodes, sLen, sPos, ePos, bNum, locDim)

    if (id==0) then   ! receive data
       if (bNum(1)==0) then   ! copy local data
           do j = 1, sN(level)
              pos = sPos(1)+(j-1)*sLen 
              pos1 = (j-1)*locDim(1)+1
              len = locDim(1)
              gaData(pos:pos+len-1)=locData(pos1:pos1+len-1)
           end do
       else 
           pos = sPos(1)
           len = ePos(1)-sPos(1)+1
           gaData(pos:pos+len-1)=locData(1:len)
       end if       

       do i=1, nNodes-1
          if (bNum(i+1)==0) then
              do j = 1, sN(level)
                 pos = sPos(i+1)+(j-1)*sLen 
                 len = locDim(i+1)
                 call MPI_Recv(gaData(pos),len,MPI_DOUBLE_PRECISION,   &
                               i, DAT_TAG, MPI_COMM_WORLD,status, ierr)
              end do
          else 
              pos = sPos(i+1)
              len = ePos(i+1)-sPos(i+1)+1
              call MPI_Recv(gaData(pos),len, MPI_DOUBLE_PRECISION,     &
                             i, DAT_TAG, MPI_COMM_WORLD, status, ierr)
          end if
       end do
    else
       if (bNum(id+1)==0) then
          do j = 1, sN(level)
             pos = (j-1)*locDim(id+1)+1
             len = locDim(id+1)
             call MPI_Send(locData(pos),len,MPI_DOUBLE_PRECISION,   &
                            0, DAT_TAG, MPI_COMM_WORLD, ierr)
          end do
       else 
           pos = 1
           len = ePos(id+1)-sPos(id+1)+1
           call MPI_Send(locData(pos),len, MPI_DOUBLE_PRECISION,     &
                         0, DAT_TAG, MPI_COMM_WORLD, status, ierr)
       end if
       
    end if
    
end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Perform direct Sending/Recving data between 2 layers:            c
!c           initial level:level1,  final level:level2                    c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine directRS(sF,sN, nNodes, id, level1, N1, dat1, level2, N2, dat2)
   implicit none
   include 'mpif.h'
   integer, intent(IN) :: sF, sN(sF), nNodes, id
   integer, intent(IN) :: level1, level2, N1, N2
   double precision, intent(IN)  :: dat1(N1)
   double precision, intent(OUT) :: dat2(N2)

   double precision :: dat0(N1)

   integer :: info, ierr, status(MPI_STATUS_SIZE)
   integer :: i, j, k, i1, i2
   integer :: sendNum,recvNum,getRecvNodesNum,getSendNodesNum, rNum1, rNum2
   integer :: bNum1(nNodes),bNum2(nNodes),locDim1(nNodes),locDim2(nNodes)
   integer(KIND=MPI_OFFSET_KIND) :: sLen1, sPos1(nNodes), ePos1(nNodes)
   integer(KIND=MPI_OFFSET_KIND) :: sLen2, sPos2(nNodes), ePos2(nNodes)
   integer, allocatable :: nInd1(:),lenInd1(:),locInd1(:),gInd1(:),req1(:), gridInd(:)
   integer, allocatable :: nInd2(:),lenInd2(:),locInd2(:),gInd2(:),req2(:)

   call getPosition(level1,sF,sN,nNodes,sLen1,sPos1,ePos1,bNum1,locDim1)
   call getPosition(level2,sF,sN,nNodes,sLen2,sPos2,ePos2,bNum2,locDim2)

   sendNum = getSendNodesNum(id,nNodes,sN(level1),sLen1,sPos1,ePos1,bNum1, &
                  sN(level2),sLen2,sPos2,ePos2,bNum2)
   recvNum = getRecvNodesNum(id,nNodes,sN(level1),sLen1,sPos1,ePos1,bNum1, &
                  sN(level2),sLen2,sPos2,ePos2,bNum2)   

   allocate(nInd1(sendNum),lenInd1(sendNum),locInd1(sendNum),gInd1(sendNum), &
            nInd2(recvNum),lenInd2(recvNum),locInd2(recvNum),gInd2(recvNum), &
            req1(sendNum), req2(recvNum),gridInd(sendNum),stat=info)
  
   call getSendNodesInd(id,nNodes,sN(level1),sLen1,sPos1,ePos1,bNum1,locDim1, &
            sN(level2),sLen2,sPos2,ePos2,bNum2,locDim2,sendNum,nInd1,lenInd1, &
            locInd1,gInd1,gridInd)   

   call getRecvNodesInd(id,nNodes,sN(level1),sLen1,sPos1,ePos1,bNum1,locDim1, &
            sN(level2),sLen2,sPos2,ePos2,bNum2,locDim2,recvNum,nInd2,lenInd2, &
            locInd2,gInd2)   
   ! receiving data
   rNum2 = 0
   do i = 1, recvNum
      if (nInd2(i)/=id) then                                     
         rNum2 = rNum2 + 1
         call MPI_IRecv(dat2(locInd2(i)),lenInd2(i),MPI_DOUBLE_PRECISION,   &
                  nInd2(i),gInd2(i),MPI_COMM_WORLD, req2(rNum2),ierr)         
      end if
   end do

   dat0(1:N1)=0.0D0

!   if ((level1==3).AND.(level2==4)) print *, 'id=',id,gridInd

   ! sending data
   rNum1 = 0
   if (bNum1(id+1)==0) then   ! grid sending
      
      do k = 1, sN(level1)
         i1 = (k-1)*locDim1(id+1)+1
         i2 = k*locDim1(id+1) 
         dat0(i1:i2)=dat1(i1:i2)      

         do i = 1, sendNum     
            if (gridInd(i)/=k) cycle
            if (nInd1(i)==id) then   ! local copy
               do j = 1, recvNum
                  if ((nInd2(j)==id).AND.(gInd1(i)==gInd2(j)))   &
                    dat2(locInd2(j):locInd2(j)+lenInd2(j)-1) =  &
                    dat0(locInd1(i):locInd1(i)+lenInd1(j)-1)
               end do
            else        ! send data
               rNum1 = rNum1 + 1
               call MPI_ISend(dat0(locInd1(i)),lenInd1(i),MPI_DOUBLE_PRECISION,   &
                  nInd1(i),gInd1(i),MPI_COMM_WORLD, req1(rNum1), ierr)
            end if
         end do

      end do
   else    ! seq sending
      dat0(1:N1)=dat1(1:N1)      
      do i = 1, sendNum     
         if (nInd1(i)==id) then   ! local copy
            do j = 1, recvNum
               if ((nInd2(j)==id).AND.(gInd1(i)==gInd2(j)))   &
                  dat2(locInd2(j):locInd2(j)+lenInd2(j)-1) =  &
                  dat0(locInd1(i):locInd1(i)+lenInd1(j)-1)
            end do
         else        ! send data
            rNum1 = rNum1 + 1
            call MPI_ISend(dat0(locInd1(i)),lenInd1(i),MPI_DOUBLE_PRECISION,   &
                  nInd1(i),gInd1(i),MPI_COMM_WORLD, req1(rNum1), ierr)
         end if
      end do
   end if


   do i = 1, rNum1
      call MPI_WAIT(req1(i), status, ierr)
   end do

   do i = 1, rNum2
      call MPI_WAIT(req2(i),status, ierr) 
   end do

   deallocate(nInd1,lenInd1,locInd1,gInd1,req1,nInd2,lenInd2,locInd2,gInd2,req2,gridInd)   
end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
