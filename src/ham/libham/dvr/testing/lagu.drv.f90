!c
!c  Testing program for lagu_podvr
!c  Functions: lagu_podvr, podvr have been tested.
!c

program lagu_podvr_test
  implicit none
  double precision, parameter :: MASS_FACTOR=1822.889D0

  integer, parameter :: N = 10 
  integer, parameter :: M = 5

  double precision, parameter :: x = 38.0D0
  double precision, parameter :: mass = 10.09D0*MASS_FACTOR
  
  double precision :: X0(M), H0(M,M), P0

  external :: pot

  integer :: lagu_podvr

  print *, 'Original Grid Points:', N
  print *, 'Final Grid Points:', M
  print *, 'Coordinate Range:[',-x, x,']'
  print *, 'Mass:', mass
  if( lagu_podvr(.FALSE., N, mass, M, POT,  X0, H0) /= 0) THEN
     print *,' error in  PODVR'
  end if

  print *, 'X node coordinates:'
  print *, X0

end 


subroutine POT(N, X, RES)
  implicit none
  double precision, parameter :: E0 = 9.0*24.743267D0  ! in cm^-1
  double precision, parameter :: A0 = 1.0D5              ! in a0

  integer, intent(IN)  :: N
  double precision, dimension(N), intent(in)  :: X
  double precision, dimension(N), intent(out) :: RES 
 
  integer :: i
  double precision :: x6

  do i = 1, N
     if (x(i) == 0.0D0) then
        RES(i) =  - E0
      else
         x6 = X(i)**6     
         RES(i) =  E0*(exp(-A0/X6) - 1)
       end if
  END DO

end

