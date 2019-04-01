!
! Testing for chevPoly, chevPolys 
!
PROGRAM TEST_Pjm
   implicit none

   integer :: J, M, N0
   double precision, allocatable :: x(:),y(:), w(:), yx(:,:)
   integer :: I, OPT = 1
   integer :: k1, k2
   double precision :: tmp

   DO WHILE (OPT /= 0)
       PRINT *, 'Input N0 (# of Gauss-Legendre integration points) '
       read  *, N0
       if ((N0 < 1) )   cycle

       allocate(x(N0), w(N0), y(N0), yx(N0, N0))
       call YjNode(N0, y)
       call YjNodes(N0, x, w)
       print *
       print *, 'Abscissas  and Weights for Gauss-Legendre Integration [-1, 1] : N=', N0
       print *, 'Abscissas:',x(1:N0)
       print *, 'Error:', X(1:N0)-Y(1:N0)
       print *, 'Weights:', w 

       if (.false.) then       
       print *, 'Test orthogonality of Associated Legendre Polynomials'
       do j = 1, N0
          print *
          print *, 'J=', j
          do m=0, j
              call YjmPolys(N0, j, m, x, yx)
              do k1 = 1, (j-m)+1
                 do k2 = 1,k1 
                  tmp = sum(yx(1:N0, k1)*yx(1:N0, k2)*w(1:N0))
                  if (k1 == k2) then
                     print *, 'M=',m, ' j1=',(m+k1-1),' j2=',(m+k2-1), 'Err :', abs(tmp-1.0D0)
                  else
                     print *, 'M=',m, ' j1=',(m+k1-1),' j2=',(m+k2-1), 'Err :', abs(tmp)
                  end if
                 end do
              end do
          end do
       end do
       end if
       deallocate (x, w, yx, y)
       PRINT *, 'INPUT OPTION for loop: O FOR EXIT'
       READ *, OPT
    END DO
END
