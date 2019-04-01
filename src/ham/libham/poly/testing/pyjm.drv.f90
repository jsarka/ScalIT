!
! Testing for chevPoly, chevPolys 
!
PROGRAM TEST_Pjm
   implicit none

   integer :: J, M, N0
   double precision :: x
   double precision, allocatable :: ch1(:),ch2(:)
   integer :: I, OPT = 1

   DO WHILE (OPT /= 0)
       PRINT *, 'Input Jmax( Pjm and Yjm )'
       read  *, J
       if (j<0)   cycle

       print *
       print *, 'Input x values: [-1, 1] : X size=', N0
       read *, x 

       N0 = (j+1)*(j+2)/2       
       allocate(ch1(N0), ch2(N0))
       call AllYjmPoly(j, x, ch1)
       call AllPjmPoly(j, x, ch2)
       print *, 'X=',x,' j=',j,' Yjm:', ch1
       print *, 'X=',x,' j=',j, ' Pjm:', ch2
       deallocate(ch1, ch2)
       
       PRINT *, 'INPUT OPTION for loop: O FOR EXIT'
       READ *, OPT
    END DO
END
