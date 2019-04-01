!ccccccccccccccccccccccccccccccccccccccccccccccccc
!c    potential for Dawes Ozone                  c
!c       Calls subroutine to convert from Jacobi c
!c       coordinates and call IMLS subroutine    c
!c       from Dawes.  Coordinates are in a.u.    c
!c       and Energy is in wavenumbers.           c
!c          - Corey Petty                        c
!ccccccccccccccccccccccccccccccccccccccccccccccccc 

double precision function potJA3(BR, lr, theta)

   use RNvar

   implicit none
   double precision, intent(IN) :: BR, lr
   DOUBLE PRECISION :: theta
   INTEGER :: numvar 
   DOUBLE PRECISION, DIMENSION(3) :: input
   INTEGER :: loadFlag = 0
   DOUBLE PRECISION, PARAMETER :: Pi = 4.*ATAN(1.)
   DOUBLE PRECISION :: rnfunc

   numneu = 50

   numw = (numinp + 1)*numneu + numneu**hidlay;
   numb = hidlay*numneu + 1;
   numvar = numw + numb;   
 
   IF (.NOT. ALLOCATED(PL)) THEN
      ALLOCATE(PL(numvar))
   END IF

   IF (loadFlag == 0) THEN
      OPEN(37, FILE='/lustre/work/corepett/ScalIT-o3well/src/systems/ozone-nn/pl.dat')
      READ(37, *) PL
      CLOSE(37)
      PRINT *, "Loaded pl.dat"
      loadFlag = 1
   END IF

   IF (theta .GT. Pi/2 .AND. theta .LE. Pi) THEN 
       theta = Pi - theta
   END IF
   IF (theta .GT. Pi .AND. theta .LE. 3*Pi/2) THEN
       theta = theta - Pi
   END IF
   IF (theta .GT. 3*Pi/2 .AND. theta .LE. 2*Pi) THEN
       theta = 2*Pi - theta
   END IF

   input(1) = lr
   input(2) = BR
   input(3) = theta

   potJA3 = rnfunc(numvar,input)

end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   1D optimized potential, using fitting parameters       c
!c        R could be in range of [0, 20au]                  c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine fitVBR(N, R0, V0)
   implicit none
   integer, intent(IN) :: N
   double precision, intent(IN)  :: R0(N)
   double precision, intent(OUT) :: V0(N)
    
   V0(1:N) = 0.0D0

end 


!cccccccccccccccccccccccccccccccccccccccccccccccc
!c      1D optimized potential for lr           c
!c     lr should be large, for example >0.5     c
!cccccccccccccccccccccccccccccccccccccccccccccccc
subroutine fitVlr(N, r0, V0)
    implicit none
    integer, intent(IN) :: N
    double precision, intent(IN)  :: r0(N)
    double precision, intent(OUT) :: V0(N)

    V0(1:N) = 0.0D0
end 

!cccccccccccccccccccccccccccccccccccccccccccccccc

