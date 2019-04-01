!
! Functions to calculate Potential of Ne3
!
subroutine potFunc3D(N, XYZ, potFun, pot)
   implicit none
   integer, intent(IN) :: N
   double precision, intent(IN) :: XYZ(3,N)
   double precision, intent(OUT) :: pot(N)
   external :: potFun

   integer :: i
   double precision :: r0(N) 

   do i = 1, N
      r0(i)=dsqrt(xyz(1,i)**2+xyz(2,i)**2+xyz(3,i)**2)
   end do
 
   call potFun(N, r0, pot)     
end 

subroutine potFunc6D(N1, XYZ1, N2, XYZ2, potFun, pot)
   implicit none
   integer, intent(IN) :: N1, N2
   double precision, intent(IN) :: XYZ1(3,N1),XYZ2(3,N2)
   external :: potFun
   double precision, intent(OUT)::pot(N1,N2)

   integer :: i, j
   double precision :: PotJA3
   double precision :: theta, r1(N1), r2(N2)

   do i = 1, N1   ! lr(i)
      r1(i) = dsqrt(XYZ1(1,i)**2+XYZ1(2,i)**2+XYZ1(3,i)**2)
   end do

   do i = 1, N2   ! BR(i)
      r2(i) = dsqrt(XYZ2(1,i)**2+XYZ2(2,i)**2+XYZ2(3,i)**2)
   end do

   do i = 1, N2
      do j = 1, N1
         theta = XYZ1(1,j)*XYZ2(1,i)+XYZ1(2,j)*XYZ2(2,i)+XYZ1(3,j)*XYZ2(3,i)        
         theta = ACOS(theta/(r1(j)*r2(i)))
         pot(j,i)=PotJA3(r2(i), r1(j), theta)
      end do
   end do

end
