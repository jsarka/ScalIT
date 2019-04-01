!ccccccccccccccccccccccccccccccccccccccccccccccccc
!c    potential for He3(LM2M2 with add-on)       c
!ccccccccccccccccccccccccccccccccccccccccccccccccc 

double precision function potJA3(BR, lr, theta)
   implicit none
   double precision, parameter  :: RMIN =4.20D0
   double precision, intent(IN) :: BR, lr, theta
    
   double precision :: x0(3)
   double precision :: cth

   cth=cos(theta);
   x0(1)=sqrt(0.25D0*lr*lr + BR*BR - BR*lr*cth);
   x0(2)=sqrt(0.25D0*lr*lr + BR*BR + BR*lr*cth);
   x0(3)=lr

   call HOO_DMBE4_PES(x0,potJA3)

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

