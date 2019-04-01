!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Subroutine written to convert the Jacobi coordinates
!  from input to Valence bond coordinates and 
!  subsequently call the Dawes Ozone PES subroutine
!
!  Written by Corey Petty 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ConvertCallPES(lr, br, gm, vpot)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: lr, br, gm 
    DOUBLE PRECISION :: r1, r2, phi
    DOUBLE PRECISION :: cgm
    DOUBLE PRECISION, DIMENSION(3) :: coords
    DOUBLE PRECISION, INTENT(OUT) :: vpot
    DOUBLE PRECISION, PARAMETER :: Pi = 3.14159265D0
    DOUBLE PRECISION, PARAMETER :: au2cm = 219474.6359029923D0  !(NIST 12/11/2013)

    cgm=COS(gm)
    r1=SQRT(0.25D0*lr*lr + br*br + lr*br*cgm) 
    r2=SQRT(0.25D0*lr*lr + br*br - lr*br*cgm)
    phi=ACOS((r1*r1 + r2*r2 - lr*lr)/(2.0D0*r1*r2))

    coords(1)=r1
    coords(2)=r2
    coords(3)=phi*180.0D0/Pi

    CALL imls( coords, vpot, 1)
  
    ! IMSL subroutine returns in cm^-1, ScalIT needs atomic units for this system. 
    vpot=vpot/au2cm

END SUBROUTINE ConvertCallPES
