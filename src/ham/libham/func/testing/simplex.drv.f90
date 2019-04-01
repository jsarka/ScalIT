!
! Testing program for minimizatation using Powell algorithm
!

program test_simplex
   implicit none
   integer, parameter :: N  = 2
   integer, parameter :: NP = 3
   double precision   :: ftol = 1.0D-6
   double precision, external :: func

   double precision, dimension(N,NP) :: P
   double precision :: fret, simplex
   
   integer :: i, iter

   P(1:N, 1:NP) = 0.0D0
   P(1,2) = 100.0D0
   P(2,3) = 100.0D0

   print *, 'Start Simplex '
   print *, 'N=', N, ' NP=', NP
   print *, 'Tolerance:', ftol
   print *, 'Initial Vertex:', P

   fret= Simplex(n, n, P, ftol, func, iter)
   print *, '# of Iteration:', iter
   print *, 'Final Vertex:', P
   print *, 'Function value:', fret

end

double precision function func(N, X)
    integer, intent(IN)  :: N
    double precision, intent(IN), dimension(N)  :: x

!    func = 2.0D0* x(1)*x(1) - 100.0D0*x(1) + 3.0D0 + &
!           x(2)*x(2) - 2.0D0*X(2) + 1.0D0


    func = 5.0D0* x(1)*x(1) - 100.0D0*x(1) + 3.0D0 + &
           x(2)*x(2) - 4.0D0*X(1)*X(2) + 1.0D0
!    print *, 'Input:',x,', Value:',func
end
