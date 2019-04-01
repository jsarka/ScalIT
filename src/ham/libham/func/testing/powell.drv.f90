!
! Testing program for minimizatation using Powell algorithm
!

program test_powell
!   implicit none
   integer, parameter :: N  = 2
   integer, parameter :: NP = 1
   double precision   :: ftol = 1.0D-6
   double precision, external :: func

   double precision, dimension(N)     :: P
   double precision, dimension(N, NP) :: Xi
   double precision :: fret
   
   integer :: i, num, powell
  
   p(1:N) = 0.0D0

   DO I = 1, NP
      Xi(1:N, I) = 0.0D0
      Xi(I,I)= 1.0D0
   END DO

   print *, 'Start POWELL '
   print *, 'N=', N, ' NP=', NP
   print *, 'Tolerance:', ftol
   print *, 'Initial Point:', P
   print *, 'Initial Direction:', Xi


   num = powell(N, NP, P, Xi, ftol, func, fret)
   print *, '# of Iteration:', num
   print *, 'Function value:', fret
   print *, 'Final Point:', P
   print *, 'Final Direction:', Xi

end

double precision function func(N, X)
    integer, intent(IN)  :: N
    double precision, intent(IN), dimension(N)  :: x

!    func = 2.0D0* x(1)*x(1) - 100.0D0*x(1) + 3.0D0 + &
!           x(2)*x(2) - 2.0D0*X(2) + 1.0D0

    func = 5.0D0* x(1)*x(1) - 100.0D0*x(1) + 3.0D0 + &
           x(2)*x(2) - 4.0D0*X(1)*X(2) + 1.0D0
end
