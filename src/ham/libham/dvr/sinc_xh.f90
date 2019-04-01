!****************************************************************
!c  Generate H matrix and DVR points for 1D sinc function       *
!c  Generating H matrix using H1_SINC1 and H1_SINC2 subroutines *
!****************************************************************
!c  General sinc.
subroutine Sinc1_XH(N, xin, mass, X0, H0)
   implicit none
   integer, intent(IN) :: N
   double precision, intent(IN) :: xin, mass
   double precision, intent(OUT):: X0(2*N+1), H0(2*N+1,2*N+1)
!c 
   call Sinc1_DVR(.false., N, xin, mass, X0, H0)
end subroutine
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
subroutine Sinc1_DVR(isCMU, N, xin, mass, X0, H0)
  logical, intent(IN)          :: isCMU
  integer, intent(IN)          :: N
  double precision, intent(IN) :: xin, mass
  double precision, intent(out):: X0(2*N+1), H0(2*N+1,2*N+1)
!c
  integer :: TOL_N, i
  double precision :: dx    
!c
  TOL_N = 2*N+1;    dx = xin/N;    X0(N+1) = 0.0D0
  do I = 1, N
     X0(N+1+I) =  dx * I
     X0(N+1-I) = -X0(N+1+I)
  end do
!c
  CALL H1_SINC1(isCMU, TOL_N, dx, mass, H0)  
!c
end subroutine
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Radial sinc.
subroutine Sinc2_XH(N, A0, B0, mass, X0, H0)
   implicit none
   integer, intent(IN) :: N
   double precision, intent(IN)  :: A0, B0, mass
   double precision, intent(OUT) :: X0(N), H0(N,N)
!c
   call Sinc2_DVR(.false., N, A0, B0, mass,X0, H0)
!c
end subroutine
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
subroutine Sinc2_DVR(isCMU, N, A0, B0, mass, X0, H0)
   implicit none
   logical, intent(IN) :: isCMU
   integer, intent(IN) :: N
   double precision, intent(IN)  :: A0, B0, mass
   double precision, intent(OUT) :: X0(N), H0(N,N)
!c
   integer :: i, i0
   double precision :: dx, xmin, xmax
!c
   xmin = MIN(A0, B0); xmax = MAX(A0, B0)
   dx = (xmax-xmin)/N; i0 = xmin/dx
   do i=1,N
      x0(i)=dx*(i0+i)
   end do
!c 
   call H1_Sinc2(isCMU, N, xmin, xmax, mass, H0)  
end subroutine
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
