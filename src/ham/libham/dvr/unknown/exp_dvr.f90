!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    DVR functions obtained from exp(imx) basis set       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine Exp_XS(Nmax, X1)
   implicit none
   integer, intent(IN) :: NMax
   double complex, intent(OUT) :: X1(2*NMax+1,2*NMax+1)

   double precision, parameter :: PI=3.1415926535897932D0
   integer :: i, j

   do i = 1, 2*NMax+1
      do j = 1, i-1
         X1(i,j)=DCMPLX(0.0D0,1.0D0/(i-j))
      end do
      X1(i,i)=PI
    end do
   
end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine Exp_XH(Nmax, mass, X1, H0)
   implicit none
   double precision, intent(IN)  :: mass 
   integer, intent(IN) :: NMax
   double complex, intent(OUT)   :: X1(2*NMax+1,2*NMax+1)
   double precision, intent(OUT) :: H0(2*NMax+1)

   double precision, parameter :: PI=3.1415926535897932D0
   integer :: i, j

   do i = 1, 2*NMax+1
      do j = 1, i-1
         X1(i,j)=DCMPLX(0.0D0,1.0D0/(i-j))
         X1(j,i)=Conjg(X1(i,j))
      end do
      X1(i,i)=PI
      H0(i)=0.5D0/mass*(i-NMax-1)**2
    end do
   
end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
