!
! Program to test CompactMat
!
program test_comp

   integer, parameter :: N = 10
   double precision :: A0(N,N), B0(N,N), C0(N,N)
   integer :: M
   double precision :: E0, H0(N),H1(N), work(3*N)
   logical :: diag

   call random_number(A0)
   C0=A0
   
   call DSYEV('N','U', N, A0, N, H0, work, 3*N, info)
   print *
   print *, ' Original Eigen value:'
   print *, H0
   print *, ' Input Cutting-Off energy:'
   read(*,*) E0
   print *, ' Ecutoff energy:',E0
   call CompactMat(E0, N, C0, M, B0)
   print *, ' Number of Compact Eigen values:', M
   if (M>0) then
      call DSYEV('N','U',M,B0,N,H1,work,3*N,info)
      print *, ' Eigen value of compact sets:'
      print *, H1(1:M)
      print *, H1(1:M)-H0(1:M)
   end if
 

end program

