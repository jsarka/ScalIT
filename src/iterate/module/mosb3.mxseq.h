!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    Subroutine to update X between layers             c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MS1X_XYZ(s1, s2, nout1, xH1, nin1, nblk1, X, nlen2, Y)
   integer, intent(IN) :: s1, s2, nin1, nout1, nblk1,nlen2
   double precision, intent(IN)    :: xH1(nout1,nout1)
   double precision, intent(INOUT) :: X(nin1,nout1,nblk1)
   double precision, intent(INOUT) :: Y(nin1,nout1,nblk1)


   double precision :: X0(nin1*nout1*nblk1),X1(nlen2)
   double precision :: Y0(nin1,nout1,nblk1),Y1(nlen2),Ytmp(nin1*nout1*nblk1)

   integer :: ierr, status(MPI_STATUS_SIZE) 
   integer :: rNum1,rNum2, ind1,ind2
   integer :: i, j, k, ii, jj

   call initLayers(s1, s2)

!ccccccccccccccccccccc   redistribute X  cccccccccccccccccccccccccccccccccc
   ! redistribute X0: receiving data
   rNum2 = 0;
   do i = 1, recvNum
      if (nInd2(i)/=id) then                                     
         rNum2 = rNum2 + 1
         call MPI_IRecv(X1(locInd2(i)),lenInd2(i),MPI_DOUBLE_PRECISION,   &
                  nInd2(i),gIndX2(i),MPI_COMM_WORLD, reqX2(rNum2),ierr)         
      end if
   end do

   ! redistribute X0: sending data
   call copyVec(nin1*nout1*nblk1,X, X0)
   rNum1 = 0
   do i = 1, sendNum     
      if (nInd1(i)==id) then   ! local copy
          do j = 1, recvNum
             if ((nInd2(j)==id).AND.(gIndX1(i)==gIndX2(j)))   &
                X1(locInd2(j):locInd2(j)+lenInd2(j)-1) =      &
                X0(locInd1(i):locInd1(i)+lenInd1(j)-1)
          end do
      else        ! send data
          rNum1 = rNum1 + 1
          call MPI_ISend(X0(locInd1(i)),lenInd1(i),MPI_DOUBLE_PRECISION,  &
                  nInd1(i),gIndX1(i),MPI_COMM_WORLD, reqX1(rNum1), ierr)
      end if
   end do


!cccccccccccccccc Y=Y+H*X, redistribute Y  cccccccccccccccccccccccccccccccc
   ! redistribute Y=H*X: receiving data
   rNum2 = 0
   do i = 1, recvNum
      if (nInd2(i)/=id) then                                     
         rNum2 = rNum2 + 1
         call MPI_IRecv(Y1(locInd2(i)),lenInd2(i),MPI_DOUBLE_PRECISION,   &
                  nInd2(i),gInd2(i),MPI_COMM_WORLD, req2(rNum2),ierr)         
      end if
   end do

   ! calculate Y=H*X, sending data   
   do k = 1, nblk1
      do i = 1, nout1
         Y0(1:nin1, i, k) = xH1(i,1)*X(1:nin1,1,k)+Y(1:nin1,i,k)
         do j=2, nout1
            Y0(1:nin1,i,k) = Y0(1:nin1,i,k)+xH1(i,j)*X(1:nin1,j, k)  
         end do
      end do
   end do

   call copyVec(nin1*nout1*nblk1,Y0,Ytmp)
   rNum1 = 0
   do ii = 1, sendNum     
      if (nInd1(ii)==id) then   ! local copy
         do jj = 1, recvNum
            if ((nInd2(jj)==id).AND.(gInd1(ii)==gInd2(jj)))    &
                 Y1(locInd2(jj):locInd2(jj)+lenInd2(jj)-1) =   &
                 Ytmp(locInd1(ii):locInd1(ii)+lenInd1(jj)-1)
         end do
      else        ! send data
         rNum1 = rNum1 + 1
         call MPI_ISend(Ytmp(locInd1(ii)),lenInd1(ii),MPI_DOUBLE_PRECISION, &
                    nInd1(ii),gInd1(ii),MPI_COMM_WORLD, req1(rNum1), ierr)
      end if     
   end do

   !cccccccccccccc Wait to finish  send/recv  cccccccccccccccccccc
   do i = 1, rNum1
      call MPI_WAIT(reqX1(i), status, ierr)
      call MPI_WAIT(req1(i),status,ierr)
   end do

   do i = 1, rNum2
      call MPI_WAIT(reqX2(i),status, ierr) 
      call MPI_WAIT(req2(i),status,ierr)
   end do
 
   call copyVec(nlen2,X1, X)
   call copyVec(nlen2,Y1, Y)

   call finalLayers()

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MS2X_DEP(s1, s2, nout1, xH1, nin1, nblk1, xDep, X, nlen2, Y)
   integer, intent(IN) :: s1, s2, nin1, nout1, nblk1,nlen2
   double precision, intent(IN)  :: xH1(nout1,nout1)
   double precision, intent(IN)  :: xDep(nin1)
   double precision, intent(INOUT) :: X(nin1,nout1,nblk1)
   double precision, intent(INOUT) :: Y(nin1,nout1,nblk1)


   double precision :: X0(nin1*nout1*nblk1),X1(nlen2)
   double precision :: Y0(nin1,nout1,nblk1),Y1(nlen2),Ytmp(nin1*nout1*nblk1)

   integer :: ierr, status(MPI_STATUS_SIZE) 
   integer :: rNum1, rNum2, ind1, ind2
   integer :: i, j, k, ii, jj

   call initLayers(s1, s2)

!ccccccccccccccccccccc   redistribute X  cccccccccccccccccccccccccccccccccc
    ! redistribute X0: receiving data
   rNum2 = 0;
   do i = 1, recvNum
      if (nInd2(i)/=id) then                                     
         rNum2 = rNum2 + 1
         call MPI_IRecv(X1(locInd2(i)),lenInd2(i),MPI_DOUBLE_PRECISION,   &
                  nInd2(i),gIndX2(i),MPI_COMM_WORLD, reqX2(rNum2),ierr)         
      end if
   end do

   ! redistribute X0: sending data
   call copyVec(nin1*nout1*nblk1,X, X0)
   rNum1 = 0
   do i = 1, sendNum
      if (nInd1(i)==id) then   ! local copy
          do j = 1, recvNum
             if ((nInd2(j)==id).AND.(gIndX1(i)==gIndX2(j)))   &
                X1(locInd2(j):locInd2(j)+lenInd2(j)-1) =      &
                X0(locInd1(i):locInd1(i)+lenInd1(j)-1)
          end do
      else        ! send data
          rNum1 = rNum1 + 1
          call MPI_ISend(X0(locInd1(i)),lenInd1(i),MPI_DOUBLE_PRECISION,  &
                  nInd1(i),gIndX1(i),MPI_COMM_WORLD, reqX1(rNum1), ierr)
      end if
   end do


!cccccccccccccccc Y=Y+H*X, redistribute Y  cccccccccccccccccccccccccccccccc   
   ! redistribute Y=H*X: receiving data
   rNum2 = 0
   do i = 1, recvNum
      if (nInd2(i)/=id) then                                     
         rNum2 = rNum2 + 1
         call MPI_IRecv(Y1(locInd2(i)),lenInd2(i),MPI_DOUBLE_PRECISION,   &
                  nInd2(i),gInd2(i),MPI_COMM_WORLD, req2(rNum2),ierr)         
      end if
   end do

   ! calculate Y=H*X, sending data   
   do k = 1, nblk1
      do i = 1, nout1
         Y0(1:nin1,i,k)=xH1(i,1)*X(1:nin1,1,k)*xDep(1:nin1)+Y(1:nin1,i,k)
         do j=2, nout1
  	    Y0(1:nin1,i,k) = Y0(1:nin1,i,k)+xH1(i,j)*X(1:nin1,j, k)*xDep(1:nin1)  
         end do
      end do
   end do

   call copyVec(nin1*nout1*nblk1,Y0,Ytmp)
  
   rNum1 = 0
   do ii = 1, sendNum     
      if (nInd1(ii)==id) then   ! local copy
         do jj = 1, recvNum
            if ((nInd2(jj)==id).AND.(gInd1(ii)==gInd2(jj)))    &
                 Y1(locInd2(jj):locInd2(jj)+lenInd2(jj)-1) =   &
                 Ytmp(locInd1(ii):locInd1(ii)+lenInd1(jj)-1)
         end do
      else        ! send data
         rNum1 = rNum1 + 1
         call MPI_ISend(Ytmp(locInd1(ii)),lenInd1(ii),MPI_DOUBLE_PRECISION,  &
                    nInd1(ii),gInd1(ii),MPI_COMM_WORLD, req1(rNum1), ierr)
      end if     
   end do


   !cccccccccccccc Wait to finish  send/recv  cccccccccccccccccccc
   do i = 1, rNum1
      call MPI_WAIT(reqX1(i), status, ierr)
      call MPI_WAIT(req1(i),status,ierr)
   end do

   do i = 1, rNum2
      call MPI_WAIT(reqX2(i),status, ierr) 
      call MPI_WAIT(req2(i),status,ierr)
   end do
 
   call copyVec(nlen2, X1, X)
   call copyVec(nlen2, Y1, Y)

   call finalLayers()

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

