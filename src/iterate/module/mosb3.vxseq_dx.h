!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    Subroutine to update X between layers             c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MS1VX_DX(s1, s2, nout1, xVi, nin1, nblk1, nlen2, X)
   integer, intent(IN) :: s1, s2, nin1, nout1, nblk1,nlen2
   double precision, intent(IN)  :: xVi(nout1,nout1,nblk1)
   double complex, intent(INOUT) :: X(nin1,nout1,nblk1)


   double complex :: X0(nlen2)
   double complex :: Y0(nin1,nout1,nblk1),Y1(nin1*nout1*nblk1)

   integer :: ierr, status(MPI_STATUS_SIZE) 
   integer :: rNum1,rNum2, ind1,ind2, lenxy
   integer :: i, j, k, ii, jj

  
   call initLayers(s1,s2)

   ! calculate Y=H*X, then redistribute Y
   ! redistribute Y=H*X: receiving data
   rNum2 = 0
   do i = 1, recvNum
      if (nInd2(i)/=id) then                                     
         rNum2 = rNum2 + 1
         call MPI_IRecv(X0(locInd2(i)),lenInd2(i),MPI_DOUBLE_COMPLEX,     &
                  nInd2(i),gInd2(i),MPI_COMM_WORLD, req2(rNum2),ierr)         
      end if
   end do

   ! calculate Y=H*X, sending data   
   do k = 1, nblk1
      do i = 1, nout1
         Y0(1:nin1, i, k) = xVi(i, 1, k) * X(1:nin1, 1, k)
         do j=2, nout1
            Y0(1:nin1,i,k) = Y0(1:nin1,i,k)+xVi(i,j,k)*X(1:nin1,j, k)  
         end do
      end do
   end do

   call copyVecCX(nin1*nout1*nblk1,Y0,Y1)
  
      ! sending data      
   rNum1 = 0
   do ii = 1, sendNum     
      if (nInd1(ii)==id) then   ! local copy
         do jj = 1, recvNum
            if ((nInd2(jj)==id).AND.(gInd1(ii)==gInd2(jj)))   &
                 X0(locInd2(jj):locInd2(jj)+lenInd2(jj)-1) =   &
                 Y1(locInd1(ii):locInd1(ii)+lenInd1(jj)-1)
         end do
      else        ! send data
         rNum1 = rNum1 + 1
         call MPI_ISend(Y1(locInd1(ii)),lenInd1(ii),MPI_DOUBLE_COMPLEX,   &
                    nInd1(ii),gInd1(ii),MPI_COMM_WORLD, req1(rNum1), ierr)
      end if     
   end do


   !cccccccccccccc Wait to finish  send/recv  cccccccccccccccccccc
   do i = 1, rNum1
      call MPI_WAIT(req1(i),status,ierr)
   end do

   do i = 1, rNum2
      call MPI_WAIT(req2(i),status,ierr)
   end do
 
   call copyVecCX(nlen2,X0, X)

   call finalLayers()

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine MS1VTX_DX(s1, s2, nout1, xVi, nin1, nblk1, nlen2, X)
   integer, intent(IN) :: s1, s2, nin1, nout1, nblk1,nlen2
   double precision, intent(IN)  :: xVi(nout1,nout1,nblk1)
   double complex, intent(INOUT) :: X(nin1,nout1,nblk1)


   double complex :: X0(nlen2)
   double complex :: Y0(nin1,nout1,nblk1),Y1(nin1*nout1*nblk1)

   integer :: ierr, status(MPI_STATUS_SIZE) 
   integer :: rNum1,rNum2, ind1,ind2, lenxy
   integer :: i, j, k, ii, jj

   call initLayers(s1,s2)

   ! calculate Y=H*X, then redistribute Y
   ! redistribute Y=H*X: receiving data
   rNum2 = 0
   do i = 1, recvNum
      if (nInd2(i)/=id) then                                     
         rNum2 = rNum2 + 1
         call MPI_IRecv(X0(locInd2(i)),lenInd2(i),MPI_DOUBLE_COMPLEX,     &
                  nInd2(i),gInd2(i),MPI_COMM_WORLD, req2(rNum2),ierr)         
      end if
   end do

   ! calculate Y=H*X, sending data   
   do k = 1, nblk1
      do i = 1, nout1
         Y0(1:nin1, i, k) = xVi(1, i, k) * X(1:nin1, 1, k)
         do j=2, nout1
            Y0(1:nin1,i,k) = Y0(1:nin1,i,k)+xVi(j,i,k)*X(1:nin1,j, k)  
         end do
      end do
   end do

   call copyVecCX(nin1*nout1*nblk1,Y0,Y1)
  
      ! sending data      
   rNum1 = 0
   do ii = 1, sendNum     
      if (nInd1(ii)==id) then   ! local copy
         do jj = 1, recvNum
            if ((nInd2(jj)==id).AND.(gInd1(ii)==gInd2(jj)))   &
                 X0(locInd2(jj):locInd2(jj)+lenInd2(jj)-1) =   &
                 Y1(locInd1(ii):locInd1(ii)+lenInd1(jj)-1)
         end do
      else        ! send data
         rNum1 = rNum1 + 1
         call MPI_ISend(Y1(locInd1(ii)),lenInd1(ii),MPI_DOUBLE_COMPLEX, &
                    nInd1(ii),gInd1(ii),MPI_COMM_WORLD, req1(rNum1), ierr)
      end if     
   end do


   !cccccccccccccc Wait to finish  send/recv  cccccccccccccccccccc
   do i = 1, rNum1
      call MPI_WAIT(req1(i),status,ierr)
   end do

   do i = 1, rNum2
      call MPI_WAIT(req2(i),status,ierr)
   end do
 
   call copyVecCX(nlen2,X0, X)

   call finalLayers()

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
