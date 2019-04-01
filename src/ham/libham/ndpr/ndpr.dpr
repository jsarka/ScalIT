!
! Get the X and H matrix using direct product-non-direct product schemeda
!
! integer function getNDPSize()
!    integer :: nR
!
!    nR = getRjSize(0)  ! # of DVR size for states=0
!    
!
! end function
!
!
! logical function getDPRXH()
!
!
! end function

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Subroutines to get X0 and H0 for the specified j        c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function getDPRH(j, N0, X1, H1)
     integer, intent(IN) :: j, N0
     double precision, intent(OUT) :: X1(N0), H1(N0, N0)

     double precision :: hmat(nMax,nMax),work(WSCALE*nMax),DH0(nMax)
     integer :: lwork, i, i1, i2, info

     lwork = WSCALE*nMax ;  getDPRH = .false.
    
     if (N0>nMax)   return

     hMat(1:nMax,1:nMax) = H0(1:nMax, 1:nMax) 
     do i=1, nMax
        hMat(i,i) = Hmat(i,i) + 0.5D0*j*(j+1)/(Mass*X0(i)*X0(i))
     end do

     call DSYEV('V', 'U', nMax, hMat, nMax, DH0, work, lwork, info)  
     if (info /= 0) return
     
     do i1 = 1, N0
        do i2 = 1, i1
	   H1(i1,i2) = SUM(hMat(1:nMax, i1)*X0(1:nMax)*hmat(1:nMax, i2))
           H1(i2,i1) = H1(i1,i2)
        end do
     end do

     call DSYEV('V', 'U', N0, H1, N0, X1, work, lwork, info)  
     if (info /= 0) return
     
     hmat(1:N0,1:N0) = h1(1:N0,1:N0)

     do i1 = 1, N0
        do i2 = 1, i1
	   H1(i1,i2) = SUM(hMat(1:N0, i1)*DH0(1:N0)*hmat(1:N0, i2))
           H1(i2,i1) = H1(i1,i2)
        end do
     end do

     getDPRH = .true.
end function
