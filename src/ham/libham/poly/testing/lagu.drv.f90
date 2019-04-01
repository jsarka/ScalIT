!
! Testing for laguPoly, laguPolys, laguNode, laguNodes 
!
PROGRAM TEST_chev
   implicit none

   integer :: M,N
   double precision, allocatable :: x0(:),x(:),w(:),lp1(:), lp2(:,:), Ln(:,:)
   integer :: I,J, OPT = 1

   DO WHILE (OPT /= 0)
       PRINT *, 'Input M (# of x value), N(Laguerre Poly Order )'
       read  *, M, N
       if ((M < 1) .or. (N<2))   cycle

       allocate(x0(M))
       print *
       print *, 'Input x values: [0, infinity] : X size=', M
       read *, x0(1:M) 

       if (M==1) then
          allocate(lp1(N+1))
          call laguPoly(N, x0, lp1)
          print *,'Laguerre Poly:x=',x0, ',N=', N,':',lp1(1:N+1) 
          
          deallocate (lp1)
       else
          allocate(lp2(M,N+1))
          call laguPolys(M,N,x0,lp2)
          DO I = 1, M               
             print *
             print *, 'X=',x0(i),'N=',N,':', lp2(i,1:N+1)
          END DO
          deallocate(lp2)
       end if

       allocate(x(2*N+1), w(2*N+1), Ln(2*N+1, N+1))
       call LaguNode(N, x)
       print *, ' N=', N
       print *, ' DVR Nodes using Direct Method:', x(1:N)
       call LaguNodes(N, x, w )
       print *
       print *, ' DVR Nodes using iterative method:',x(1:N)
       print *, ' DVR weights:', w(1:N)
       print *  

       call laguNodes(2*N+1, x, w)
       call laguPolys(2*N+1, N, x, Ln )
       print *, 'Testing orthogonality of Ln'
       do i = 1, N
          do j = 1, N
                print *, '<',i,'|',j,'>=',sum(w(1:2*N+1)*Ln(1:2*N+1,i)*Ln(1:2*N+1,j))
          end do
       end do
       deallocate (x0, x, w, ln)
       PRINT *, 'INPUT OPTION for loop: O FOR EXIT'
       READ *, OPT
    END DO
end
