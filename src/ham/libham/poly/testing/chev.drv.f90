!
! Testing for chevPoly, chevPolys 
!
PROGRAM TEST_chev
   implicit none

   integer :: M,N
   double precision, allocatable :: x(:),ch1(:),ch2(:,:),ch0(:)
   integer :: I,J, OPT = 1

   DO WHILE (OPT /= 0)
       PRINT *, 'Input M (# of x value), N(Chev Poly Order )'
       read  *, M, N
       if ((M < 1) .or. (N<2))   cycle

       allocate(x(M))
       print *
       print *, 'Input x values: [-1, 1] : X size=', M
       read *, x(1:M) 

       if (M==1) then
          allocate(ch1(N+1), ch0(N+1))
          call chevPoly(N, x, ch1)
          do i = 0, N
             ch0(i+1) = cos(i*acos(x(1)))
          end do
          do i = 0, N
              print *,'Tm:x=',x, ',N=', i,':',ch1(i+1),' DE=',ch0(i+1)-ch1(i+1) 
          end do
          deallocate (ch1, ch0)
       else
          allocate(ch2(M,N+1), ch0(N+1))
          call chevPolys(M,N,x,ch2)
          DO I = 1, M               
             do j = 0, N
                ch0(j+1) = cos(j*acos(x(i)))
             end do
             print *, 'X=',x(i),'N=',N,':', ch2(i,1:N+1)
             print *, 'Error', ch0(1:N+1)-ch2(i,1:N+1)
          END DO
          deallocate(ch2, ch0)
       end if
       deallocate (x)
       PRINT *, 'INPUT OPTION for loop: O FOR EXIT'
       READ *, OPT
    END DO
END
