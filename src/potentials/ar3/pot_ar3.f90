!
!  potential for trimeter of rare gas
!  using L-J potential
!

double precision function potJA3(br, lr, theta)
   implicit none
   double precision, parameter  :: RMIN =0.50D0
   double precision, intent(IN) :: BR, lr, theta
    
   double precision :: lx, BRx, BRy
   double precision :: r1_2, r2_2, r3_2, r1_6, r2_6, r3_6

   lx = 0.5D0*lr
   BRx = BR*cos(theta); BRy = BR*sin(theta);
   r1_2  = lr * lr;
   r2_2  = (BRx-lx)*(BRx-lx) + BRy*BRy
   r3_2  = (BRx+lx)*(BRx+lx) + BRy*BRy

   if (r1_2<RMIN) r1_2 = RMIN
   if (r2_2<RMIN) r2_2 = RMIN
   if (r3_2<RMIN) r3_2 = RMIN

   r1_2 = 1.0D0/r1_2; r2_2 = 1.0D0/r2_2; r3_2 = 1.0D0/r3_2;
   r1_6 = r1_2*r1_2*r1_2;  r2_6 = r2_2*r2_2*r2_2
   r3_6 = r3_2*r3_2*r3_2;

   potJA3 = (r1_6*r1_6 - r1_6) + (r2_6*r2_6 - r2_6) + (r3_6*r3_6 - r3_6)
   potJA3 = 4.0D0*potJA3

end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   1D optimized potential, using fitting parameters       c
!c        R could be in range of [0, infinity)              c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine fitVBR(N, R0, V0)
   implicit none
   double precision, parameter:: Coeff(32)                                    &
      = (/ -25.76770842336366D0,  44.42081721973262D0, -50.01114123911642D0,  & 
            -1.96307362468141D0, -11.36877790896871D0,   0.93805414891994D0,  &
            -1.75669014587530D0,   1.88353784555610D0,   3.11511139000352D0,  &
             0.83089765422437D0,  -1.03913967831511D0,  -3.38116671332510D0,  &
            -0.34546710402757D0,   1.70265867428934D0,  -0.55194129629639D0,  &
             0.25344030367140D0,  -1.08381086076204D0,  19.91713085284756D0,  &
             8.09775638753977D0,   9.27546129753706D0,  12.27595577031276D0,  &
            -1.29397596036303D0,   1.93829501363316D0,   1.74605203237829D0,  &
           -13.77747239433938D0,   1.74295584871104D0,  -0.38533510994948D0,  &
             0.09160659533092D0,  -0.51890338976101D0,   0.43172233025939D0,  &
             0.40071376266585D0,  -0.74021862992175D0 /)
    integer, intent(IN) :: N
    double precision, intent(IN)  :: R0(N)
    double precision, intent(OUT) :: V0(N)
    
    double precision :: X1(N), X2(N), X3(N)

    X1(1:N) = exp(coeff(1)*R0(1:N)*R0(1:N)+coeff(2)*R0(1:N)+coeff(3)) 
    call polyval(4, coeff(13:17), N, R0, X2);
    call mystep(coeff(4:7), N, R0, V0)
    V0(1:N) = V0(1:N)*(X1(1:N)+X2(1:N));
 
    X1(1:N) = exp((coeff(18)*(R0(1:N)+coeff(19)))**(-6));    
    call mystep(coeff(8:11), N, R0, X2)
    V0(1:N) = V0(1:N) + X1(1:N)*X2(1:N); 

    x1(1:N) = exp((coeff(20)*(R0(1:N)+coeff(21)))**(-6)); 
    call polyVal(4, coeff(28:32), N, R0, X2)
    X1(1:N) = X1(1:N) + X2(1:N)
    call mymidFilter(coeff(22:27), N, R0, X2);

    V0(1:N) = V0(1:N) + X1(1:N)*X2(1:N) + coeff(12);   
end 


!cccccccccccccccccccccccccccccccccccccccccccccccc
!c      1D optimized potential for lr           c
!c     lr should be large, for example >0.5     c
!cccccccccccccccccccccccccccccccccccccccccccccccc
subroutine fitVlr(N, r0, V0)
    implicit none
    double precision, parameter :: coeff(9) =(/                           &
         -48.77170394213497D0, -6.67558890379387D0,  2.28453344352343D0,  & 
          52.81368877933721D0,  0.52342807916784D0,  4.57934885619104D0,  &
           2.35031076831229D0, -0.59408702095415D0, -0.88798240321165D0  /)
    integer, intent(IN) :: N
    double precision, intent(IN)  :: r0(N)
    double precision, intent(OUT) :: V0(N)

    double precision :: X1(N), X2(N), r1(N)
    integer :: i
    double precision,parameter :: rmin=0.5D0
    
    r1(1:N) = r0(1:N)
    do i=1,N
       if (r0(i)<rmin) r1(i)=rmin
    end do

    X1(1:N) = r1(1:N)**(-12) - r1(1:N)**(-6);
    call mystep(coeff(1:4), N, r1, X2)
    V0 = X1(1:N)*X2(1:N) ;
    
    X1(1:N) = exp((10.0D0*r1(1:N))**(-6));
    call mystep(coeff(5:8), N, r1, X2)
    
    V0(1:N)=V0(1:N) + X1(1:N)*X2(1:N) + coeff(9);

end 

!cccccccccccccccccccccccccccccccccccccccccccccccc
!c      1D optimized potential for lr           c
!c     lr should be large, for example >0.5     c
!cccccccccccccccccccccccccccccccccccccccccccccccc
subroutine fitVlr0(N, r0, V0)
    implicit none
    double precision, parameter :: coeff(9) =(/                           &
         -48.77170394213497D0, -6.67558890379387D0,  2.28453344352343D0,  &
          52.81368877933721D0,  0.52342807916784D0,  4.57934885619104D0,  &
           2.35031076831229D0, -0.59408702095415D0, -0.88798240321165D0  /)
    integer, intent(IN) :: N
    double precision, intent(IN)  :: r0(N)
    double precision, intent(OUT) :: V0(N)

    double precision :: X1(N), X2(N)

    X1(1:N) = r0(1:N)**(-12) - r0(1:N)**(-6);
    call mystep(coeff(1:4), N, r0, X2)
    V0 = X1(1:N)*X2(1:N) ;

    X1(1:N) = exp((10.0D0*r0(1:N))**(-6));
    call mystep(coeff(5:8), N, r0, X2)

    V0(1:N)=V0(1:N) + X1(1:N)*X2(1:N) + coeff(9);

end

