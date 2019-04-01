!ccccccccccccccccccccccccccccccccccccccccccccccccc
!c    potential for SO2 -Corey Petty             c
!ccccccccccccccccccccccccccccccccccccccccccccccccc 

double precision function potJA3(BR, lr, theta)
   implicit none
   double precision, intent(IN) :: BR, lr, theta
   double precision :: potential
   double precision , dimension(2) :: x0
   double precision :: phi
   double precision :: cth
   double precision :: sth
   double precision :: Pi
   character :: inp*3, mod*2

   !inp="DTT"    !different data sets, see PES
   inp="frr"
!!!!!!convert from jacobi to valance bond coordinates!!!!!!!!!!!!
   Pi=3.14159265
   cth=cos(theta)
   sth=sin(theta)
   x0(1)=sqrt(0.25D0*lr*lr + BR*BR + lr*BR*cth)
   x0(2)=sqrt(0.25D0*lr*lr + BR*BR - lr*BR*cth)
   phi=asin(0.5D0*lr*sth/x0(2)) + asin(0.5D0*lr*sin(Pi-theta)/x0(1))

!!!!!!convert from a.u. to angstroms!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   x0(1)=x0(1)*0.529177249D0
   x0(2)=x0(2)*0.529177249D0
   

!!!!!!call potential energy surface subroutine!!!!!!!!!!!!!!!!!!!
!!!!!!corresponds to potff3.f file!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if(phi < 0.0D0) print*, "!! THE ANGLE IS NEGATIVE!!!"
   if(phi> Pi)    print*, "!! THE ANGLE IS > PI !!!!!!"
   if(theta .GE. 1.0D0 .AND. theta .LE. 2.0D0) then
      call potff3(x0,phi,potential,inp,mod)
      potential=potential/27.2113961D0
   else    
      potential=3.0D0/27.2113961D0
   end if  
!!!!!!convert back to a.u. energies from eV!!!!!!!!!!!!!!!!!!!!!!
   if(potential <= 3.0D0/27.2113961D0 .AND. potential > 0.0D0) then
         potJA3=potential
   else 
         potJA3=3.0D0/27.2113961D0
   end if
!   potJA3=0.0D0
!   print*, "The Potential is: " , potJA3
   RETURN
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

