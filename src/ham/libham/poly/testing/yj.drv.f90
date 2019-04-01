!
! Testing for chevPoly, chevPolys 
!
PROGRAM TEST_Yj
   implicit none

   integer :: M,N
   double precision, allocatable :: x(:),ch1(:),ch2(:,:)
   integer :: I, OPT = 1

   DO WHILE (OPT /= 0)
       PRINT *, 'Input M (# of x value), j(Norm. Legendre Poly Order )'
       read  *, M, N
       if ((M < 1) .or. (N<2))   cycle

       allocate(x(M))
       print *
       print *, 'Input x values: [-1, 1] : X size=', M
       read *, x(1:M) 

       if (M==1) then
          allocate(ch1(N+1))
          call YjPoly(N, x, ch1)
          print *, 'Yj: x=',x, ', Y(j=', N,'):',ch1(:)
          deallocate (ch1)
       else
          allocate(ch2(M,N+1))
          call YjPolys(M,N,x,ch2)
          DO I = 1, M        
             print *, 'X=',x(i),'Y(j=',N,'):', ch2(i,:)
          END DO
          deallocate(ch2)
       end if
       deallocate (x)
       PRINT *, 'INPUT OPTION for loop: O FOR EXIT'
       READ *, OPT
    END DO
END
