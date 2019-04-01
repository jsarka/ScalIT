!
! Testing for chevPoly, chevPolys 
!
PROGRAM TEST_Pjm
   implicit none

   integer :: J, M, N0
   double precision, allocatable :: x(:),ch1(:),ch2(:,:)
   integer :: I, OPT = 1

   DO WHILE (OPT /= 0)
       PRINT *, 'Input N0 (# of x value), M, J( Associated Legendre Pjm )'
       read  *, N0, M, J
       if ((N0 < 1) .or. (j<m))   cycle

       allocate(x(N0))
       print *
       print *, 'Input x values: [-1, 1] : X size=', N0
       read *, x(1:N0) 

       if (N0==1) then
          allocate(ch1(J))
          call PjmPoly(j, m, x, ch1)
          print *, 'Pj: x=',x, ', (j,m)=', j,M,'):',ch1(1:j-m+1)
          deallocate (ch1)
       else
          allocate(ch2(N0,J))
          call PjmPolys(N0,j, m, x, ch2)
          DO I = 1, N0       
             print *, 'X=',x(i),'(j,m)=',j,m,':', ch2(i,1:(j-m)+1)
          END DO
          deallocate(ch2)
       end if
       deallocate (x)
       PRINT *, 'INPUT OPTION for loop: O FOR EXIT'
       READ *, OPT
    END DO
END
