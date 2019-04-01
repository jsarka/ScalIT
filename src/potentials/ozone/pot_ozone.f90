!ccccccccccccccccccccccccccccccccccccccccccccccccc
!c    potential for Dawes Ozone                  c
!c       Calls subroutine to convert from Jacobi c
!c       coordinates and call IMLS subroutine    c
!c       from Dawes.  Coordinates are in a.u.    c
!c       and Energy is in wavenumbers.           c
!c          - Corey Petty                        c
!ccccccccccccccccccccccccccccccccccccccccccccccccc 

double precision function potJA3(BR, lr, theta)
   implicit none
   double precision, intent(IN) :: BR, lr, theta
    
   call ConvertCallPES(lr, BR, theta, potJA3)

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

