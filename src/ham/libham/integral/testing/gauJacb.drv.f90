!
! Testing for DVR points using Lagauerr polynomials 
!  NOTE: This seems to be just a copy of lagauerrs test driver
!     and wont actually test anything.  Should give a guide on
!     how to make it work though.  -C.P.
!
PROGRAM TEST_Gau_Lagaueer
   implicit none

   integer :: N0
   double precision, allocatable :: x0(:),x1(:),x2(:),w1(:),w2(:)

       PRINT *, 'Input N0 (# of Gauss-Laguerre integration points) '
       read  *, N0

       allocate(x0(N0), x1(N0), x2(N0), w1(N0), w2(N0))
       call LaguX0(N0, x0)

       call LaguX0W0(N0, x1, w1)

       call GauLaguX0W0(N0, X2, W2)

       print *
       print *, 'Abscissas and Weights for Gauss-Laguerre Integration [-1, 1]:N=', N0
       print *, 'Abscissas:'
       print *, x1(1:N0)
       print *, 'Weights:'
       print *, w1 

       print *
       print *, 'X0 error:,LaguX0-LaguX0W0'
       print *, X0-X1

       print *
       print *, 'Error: LaguX0-GauLaguX0'
       print *,  X0-X2

       print *
       print *, ' W0 Error:'
       print *, w1-w2

       deallocate (x0, x1, x2, w1, w2)

END
