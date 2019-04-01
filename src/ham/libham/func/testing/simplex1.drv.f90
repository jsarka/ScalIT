!
! Testing program for minimizatation using Powell algorithm
!

program test_ameoba
   implicit none
   integer, parameter :: N  = 2
   integer, parameter :: NP = 3
   double precision   :: ftol = 1.0D-6
   double precision, external :: func

   double precision, dimension(NP,N) :: P
   double precision,dimension(NP) :: xi
   double precision :: fret
   
   integer :: i, iter

   P(1:NP, 1:N) = 0.0D0
   P(2,1) = 100.0D0
   P(3,2) = 100.0D0

   xi(1) = func(p(1,1:n))
   xi(2) = func(p(2,1:n))
   xi(3) = func(p(3,1:n))
   print *, 'Start Simplex '
   print *, 'N=', N, ' NP=', NP
   print *, 'Tolerance:', ftol
   print *, 'Initial Vertex:', P
   print *, 'Initial value:', xi

   call amoeba(P, xi, np, n, n, ftol, func, iter)
   print *, '# of Iteration:', iter
   print *, 'Final Vertex:', P
   print *, 'Function value:', xi

end

double precision function func(X)
!    integer, intent(IN)  :: N
    double precision, intent(IN), dimension(2)  :: x

!    func = 2.0D0* x(1)*x(1) - 100.0D0*x(1) + 3.0D0 + &
!           x(2)*x(2) - 2.0D0*X(2) + 1.0D0


    func = 5.0D0* x(1)*x(1) - 100.0D0*x(1) + 3.0D0 + &
           x(2)*x(2) - 4.0D0*X(1)*X(2) + 1.0D0
!    print *, 'Input:',x,', Value:',func
end
