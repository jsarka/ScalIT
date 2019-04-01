!
! Testing for chevPoly, chevPolys 
!
PROGRAM TEST_Pj
   implicit none

   integer :: M,N
   double precision, allocatable :: x(:),ch1(:),ch2(:,:)
   integer :: I, OPT = 1

   DO WHILE (OPT /= 0)
       PRINT *, 'Input M (# of x value), L(Legendre Poly Order )'
       read  *, M, N
       if ((M < 1) .or. (N<2))   cycle

       allocate(x(M))
       print *
       print *, 'Input x values: [-1, 1] : X size=', M
       read *, x(1:M) 

       if (M==1) then
          allocate(ch1(N+1))
          call PjPoly(N, x, ch1)
          print *, 'Pj: x=',x, ', P(l=', N,'):',ch1(:)
          deallocate (ch1)
       else
          allocate(ch2(M,N+1))
          call PjPolys(M,N,x,ch2)
          DO I = 1, M        
             print *, 'X=',x(i),'P(l=',N,'):', ch2(i,:)
          END DO
          deallocate(ch2)
       end if
       deallocate (x)
       PRINT *, 'INPUT OPTION for loop: O FOR EXIT'
       READ *, OPT
    END DO
END
