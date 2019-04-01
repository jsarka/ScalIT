! Testing exp_drv, dvrcx0, dvrcx0SF, dvrcx0RF
program test_dvrcx
   implicit none
   integer, parameter :: NMAX=3, N=2*NMAX+1, M=7
   double precision :: mass = 1.0D0   
   external :: pot
   double complex :: X1(N,N), H1(N), H0(M,M)
   double precision :: X0(M),E0(M)

   call Exp_XH(NMAX, mass, X1, H1)

   call DVRCX0(N,X1,H1,pot,M,X0,E0,H0)

   print *, ' DVR points:', X0
   print *, ' Eigen values:',E0
   print *, ' H matrix:',H0

end


subroutine pot(N,X,pv)
   integer :: N
   double precision :: X(N), pv(N)   

   pv(1:N) = 0.5*X(1:N)*X(1:N)
end


