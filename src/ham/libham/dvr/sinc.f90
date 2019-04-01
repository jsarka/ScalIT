!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Define some functions to calculate H1 matrix                            c
!c  There are functions to calculate H1 matrix using SINC functions         c
!c  H1 FOR SINC FUNCTION: See JCP 96(3), 1982 (1992).                       c
!c  It doesn't work for angular variables                                   c
!c  The energy is in cm^-1 when isCMU is .true., otherwise in Hartree unit  c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c  The DVR or grid point representation of the kinetic energy
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         H1_SINC1: FOR [-INFINITY, INFINITY] interval            c
!c                                 PI^2 / 3                        c
!c  H1=1/(2*mass*dx^2)*(-1)^(I-I')                                 c
!c                                 2/(I-I')**2                     c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
  subroutine H1_SINC1(isCMU, N, dx, mass, H1)
      implicit none
      logical, intent(IN)          :: isCMU
      integer, intent(IN)          :: N
      double precision, intent(IN) :: dx, mass
      double precision, intent(OUT):: H1(N,N)
!c
      double precision, parameter :: PI = 3.14159265358979323846D0
      double precision, parameter :: PI3 = PI**2/3.0D0
      double precision, parameter :: HARTREE_CM  = 8.06554477D3 * 27.2114D0
!c      
      integer :: I, J, IJ
      double precision :: SCALE
!c
      if (isCMU) then
          SCALE = 0.50D0/(dx*dx*mass) * Hartree_CM   !..in cm inverse
      else 
          SCALE = 0.50D0/(dx*dx*mass)                !..in Hartree
      end if
!c
      do I = 1, N
         do J = 1, I-1
            IJ = I - J
            if (IJ/2 * 2 == IJ) then
               H1(I, J) = 2.0D0 / (IJ*IJ)
            else
               H1(I, J) = -2.0D0 / (IJ*IJ)
            end if
            H1(J, I) = H1(I, J)
         end do
      end do
!c      
      do I = 1, N
         H1(I,I) = PI3
      end do
!c
      H1(1:N, 1:N) = SCALE * H1(1:N, 1:N)
!c
  end subroutine
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  This is appropriate for a radial coordinate r                  c
!c    H1_SINC2: FOR [0, INFINITY] interval                        c
!c                           PI^2 / 3  - 1/(2*I**2)                c
!c H1 = 1/(2*M*dx^2)*(-1)^(I-I')                                   c
!c                           2/(I-I')**2 - 2/(I+I')**2             c
!c  A0: Initial grid point; B0: End grid ponit                     c 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
  subroutine H1_SINC2(isCMU, N, A0, B0, mass, H1)
      implicit none
      logical, intent(IN) :: isCMU
      integer, intent(IN) :: N
      double precision, intent(IN) :: A0, B0, mass
      double precision, intent(OUT) :: H1(N,N)
!c
      double precision, parameter :: PI = 3.14159265358979323846D0
      double precision, parameter :: PI3 = PI**2/3.0D0
      double precision, parameter :: HARTREE_CM  = 8.06554477D3 * 27.2114D0
!c      
      integer :: I, J, I0
      integer :: IJ1, IJ2
      double precision :: dx, SCALE
!c      
      dx = (B0-A0)/N
      I0 = A0 / dx
!c
      if (isCMU) then
         SCALE = 0.50D0/(dx*dx*mass)*Hartree_CM    !..in cm inverse
      else
         SCALE = 0.50D0/(dx*dx*mass)               !..in Hartree
      end if
!c
      do I = 1, N
         do J = 1, I-1
            IJ1 = I - J
            IJ2 = I + J + I0 + I0
            H1(I, J) = 2.0D0/(IJ1*IJ1) - 2.0D0/(IJ2*IJ2)
            if ((IJ1/2 * 2) /=  IJ1) then            
               H1 (I, J) = - H1(I, J)
            end if
            H1(J, I) = H1(I, J)
         end do
      end do
!c      
      do I = 1, N
         H1(I,I) = PI3 - 0.5D0/((I+I0)*(I+I0))
      end do
!c
      H1(1:N, 1:N) = SCALE * H1(1:N, 1:N)
!c
  end subroutine
!c
!***********************************************************************
