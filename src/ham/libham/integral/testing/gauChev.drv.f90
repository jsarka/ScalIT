!
! Testing for DVR points using Chebeshev polynomials 
!
PROGRAM TEST_Gau_Chebeshev
   implicit none

   integer :: N0
   double precision, allocatable :: x0(:),x1(:),x2(:),w1(:),w2(:)

       PRINT *, 'Input N0 (# of Gauss-Legendre integration points) '
       read  *, N0

       allocate(x0(N0), x1(N0), x2(N0), w1(N0), w2(N0))
       call ChevX0(N0, x0)

       call ChevX0W0(N0, x1, w1)

       call GauChevX0W0(N0, X2, W2)

       print *
       print *, 'Abscissas and Weights for Gauss-Chebeshev Integration [-1, 1]:N=', N0
       print *, 'Abscissas:'
       print *, x1(1:N0)
       print *, 'Weights:'
       print *, w1 

       print *
       print *, 'X0 error:,ChevX0-ChevX0W0'
       print *, X0-X1

       print *
       print *, 'Error: ChevX0-GauCheveX0'
       print *,  X0+X2

       print *
       print *, ' W0 Error:'
       print *, w1-w2

       deallocate (x0, x1, x2, w1, w2)

END
