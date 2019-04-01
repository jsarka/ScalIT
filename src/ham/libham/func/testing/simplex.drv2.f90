!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                    !
!     A simple test program                                          !
!                                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM tminim
! Use minim to maximize the objective function:
! 1/{1 + (x-y)^2} + sin(pi.y.z/2) + exp[-{(x+z)/y - 2}^2]
! with respect to x, y and z.
! We will actually minimize its negative.
! The maximum occurs at x = y = z = +/-sqrt(4n+1) for any integer n.
IMPLICIT NONE
DOUBLE PRECISION :: object, p(3), simp, step(3), stopcr, var(3)
INTEGER          :: ier, iprint, iquad, maxf, nloop, nop
LOGICAL          :: first
! remember to declare the user's subroutine as EXTERNAL.
EXTERNAL objfun
! Set up starting values & step sizes.
! Parameters p(1), p(2), p(3) are x, y & z.
p(1) = 0.d0
p(2) = 1.d0
p(3) = 2.d0
step = 0.4D0
nop = 3
! Set max. no. of function evaluations = 250, print every 10.
maxf = 250
iprint = 10
! Set value for stopping criterion.   Stopping occurs when the
! standard deviation of the values of the objective function at
! the points of the current simplex < stopcr.
stopcr = 1.d-04
nloop = 6
! Fit a quadratic surface to be sure a minimum has been found.
iquad = 1
! As function value is being evaluated in double precision, it
! should be accurate to about 15 decimals.   If we set simp = 1.d-6,
! we should get about 9 dec. digits accuracy in fitting the surface.
simp = 1.d-6
! Now call MINIM to do the work.
first = .true.

DO
  CALL minim(p, step, nop, object, maxf, iprint, stopcr, nloop,   &
              iquad, simp, var, objfun, ier)
! If ier > 0, try a few more function evaluations.
  IF (ier .EQ. 0) EXIT
  IF (.NOT. first) STOP
  first = .false.
  maxf = 100
END DO
! Successful termination.
WRITE(*, 900) object, p
900 FORMAT(' Success !'/' Objective function = ', f12.6/ ' at: ', 3f12.6)
END PROGRAM tminim

SUBROUTINE objfun(p, func)
! This is the subroutine which the user must write.
! Remember that we are minimizing the negative of the function we
! really want to maximize.
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN)  :: p(3)
DOUBLE PRECISION, INTENT(OUT) :: func
!     Local variables
DOUBLE PRECISION :: half = 0.5D0, one = 1.D0, pi = 3.14159265357989D0,   &
                    two = 2.D0, x, y, z
x = p(1)
y = p(2)
z = p(3)
func = -one/(one + (x-y)**2) - SIN(half*pi*y*z) - EXP(-((x+z)/y - two)**2)
RETURN
END SUBROUTINE objfun


