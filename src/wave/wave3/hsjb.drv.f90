!
! Test hs2jb() and jb2hs() subroutines
!
program hsjb
   implicit none
   integer, parameter :: NUM=20
   double precision, parameter :: R0MAX=10.0, r1Max=5.0
   double precision, dimension(NUM) :: R0, r1, phi
   double precision :: R0_1, r1_1, phi_1, cphi, pi
   double precision :: rho, theta, chi, ththeta, c2chi, s2chi
   integer :: i
 
   call random_number(R0)
   call random_number(r1)
   call random_number(phi)
   R0=R0MAX*R0;  r1=r1*r1Max
   pi = dcos(-1.0D0)

   print *
   print *, ' Transform between Jacobi-Hyperspherical coordinates! '

   do i = 1, NUM
      print *
      print *, ' I=', i
      write(*,10) R0(i),r1(i),phi(i)
      call jb2hs(R0(i),r1(i), phi(i), rho, theta, chi)

      write(*,30) rho, theta, chi
      call hs2jb(rho,theta,chi, R0_1,r1_1, cphi)
      phi_1=acos(cphi)
      
      write(*,10) R0_1, r1_1, phi_1
      write(*,20) R0(i)-R0_1,r1(i)-r1_1,phi(i)-phi_1
   end do
    
   print *, ' Finish the testing!'
   print *

   10 format(' Jacobi:  R0=',F15.9, 2x,'  r1=  ', F15.9, 2x, ' Phi=', F15.9)
   20 format(' Diff:    R0=',F15.9, 2x,'  r1=  ', F15.9, 2x, ' Phi=', F15.9)
   30 format(' HyperSP:rho=',F15.9, 2x,' Theta=', F15.9, 2x, ' Chi=', F15.9)

end

