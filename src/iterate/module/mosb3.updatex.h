!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    Subroutine to update X between layers             c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine updateX(s1, s2, dat1, dat2)
   integer, intent(IN) :: s1, s2
   double precision, intent(IN)  :: dat1(plen(s1))
   double precision, intent(OUT) :: dat2(plen(s2))

   integer :: ierr, status(MPI_STATUS_SIZE) 
   integer :: i, j, rNum1, rNum2

   call initLayers(s1, s2)

!   if (id==rootID) then
!     print *, 'node:',recvNum, sendNum, s1, s2
!     print *, 'sLen1:',sLen(s1), ' sLen2=',sLen(s2)
!     print *, 'sPos1:', sPos(:,s1)
!     print *, 'ePos1:', ePos(:,s1)
!     print *, 'bNum1:',bNum(:,s1)
!     print *, 'locDim1:',locDim(:,s1)
!     print *, 'gInd1:',gInd1(:)
!     print *, 'gridInd:',gridInd(:)
!     print *, 'sPos2:', sPos(:,s2)
!     print *, 'ePos2:', ePos(:,s2)
!     print *, 'bNum2:',bNum(:,s2)
!     print *, 'locDim2:',locDim(:,s2)
!     print *, 'gInd2:',gInd2(:)
!
!   end if


   ! sending data
   rNum1 = 0
   do i = 1, sendNum     
      if (nInd1(i)==id) then   ! local copy
          do j = 1, recvNum
             if ((nInd2(j)==id).AND.(gInd1(i)==gInd2(j)))   &
                dat2(locInd2(j):locInd2(j)+lenInd2(j)-1) =  &
                dat1(locInd1(i):locInd1(i)+lenInd1(j)-1)
          end do
      else        ! send data
          rNum1 = rNum1 + 1
          call MPI_ISend(dat1(locInd1(i)),lenInd1(i),MPI_DOUBLE_PRECISION,   &
                  nInd1(i),gInd1(i),MPI_COMM_WORLD, req1(rNum1), ierr)
      end if
   end do

   ! receiving data
   rNum2 = 0
   do i = 1, recvNum
      if (nInd2(i)/=id) then                                     
         rNum2 = rNum2 + 1
         call MPI_IRecv(dat2(locInd2(i)),lenInd2(i),MPI_DOUBLE_PRECISION,   &
                  nInd2(i),gInd2(i),MPI_COMM_WORLD, req2(rNum2),ierr)         
      end if
   end do

   do i = 1, rNum1
      call MPI_WAIT(req1(i), status, ierr)
   end do

   do i = 1, rNum2
      call MPI_WAIT(req2(i),status, ierr) 
   end do

   call finalLayers()

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine updateX_CX(s1, s2, dat1, dat2)
   integer, intent(IN) :: s1, s2
   double complex, intent(IN)  :: dat1(plen(s1))
   double complex, intent(OUT) :: dat2(plen(s2))

   integer :: ierr, status(MPI_STATUS_SIZE)
   integer :: i, j, rNum1, rNum2

   integer,allocatable :: nInd1(:),lenInd1(:),locInd1(:),gInd1(:),req1(:),gridInd(:)
   integer,allocatable :: nInd2(:),lenInd2(:),locInd2(:),gInd2(:),req2(:)

   call initLayers(s1, s2)

   ! sending data
   rNum1 = 0
   do i = 1, sendNum     
      if (nInd1(i)==id) then   ! local copy
          do j = 1, recvNum
             if ((nInd2(j)==id).AND.(gInd1(i)==gInd2(j)))   &
                dat2(locInd2(j):locInd2(j)+lenInd2(j)-1) =  &
                dat1(locInd1(i):locInd1(i)+lenInd1(j)-1)
          end do
      else        ! send data
          rNum1 = rNum1 + 1
          call MPI_ISend(dat1(locInd1(i)),lenInd1(i),MPI_DOUBLE_COMPLEX,   &
                  nInd1(i),gInd1(i),MPI_COMM_WORLD, req1(rNum1), ierr)
      end if
   end do

   ! receiving data
   rNum2 = 0
   do i = 1, recvNum
      if (nInd2(i)/=id) then                                     
         rNum2 = rNum2 + 1
         call MPI_IRecv(dat2(locInd2(i)),lenInd2(i),MPI_DOUBLE_COMPLEX,   &
                  nInd2(i),gInd2(i),MPI_COMM_WORLD, req2(rNum2),ierr)         
      end if
   end do

   do i = 1, rNum1
      call MPI_WAIT(req1(i), status, ierr)
   end do

   do i = 1, rNum2
      call MPI_WAIT(req2(i),status, ierr) 
   end do

   call finalLayers()

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
