!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    Calculate full potential in Jacobi Coordinates     c
!c        (BR, lr1, lr2, theta1, theta2, phi)            c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! lr1=rCC, lr2=rHH, BR=rCC-HH
!
subroutine potJA4(BR, lr1, lr2, N1, N2, N3,th1,th2,phi,V) 
   implicit none
   double precision, intent(IN)  :: BR, lr1, lr2
   integer, intent(IN) :: N1, N2, N3
   double precision, intent(IN)  :: th1(N1),th2(N2), phi(N3)
   double precision, intent(OUT) :: V(N3, N1, N2)

   double precision :: cth1(N1),sth1(N1),cth2(N2),sth2(N2),cphi(N3),sphi(N3)
   double precision :: r(6)
   integer :: i1, i2, i3

   cth1(1:N1)=cos(th1(1:N1)); sth1(1:N1)=sin(th1(1:N1))
   cth2(1:N2)=cos(th2(1:N2)); sth2(1:N2)=sin(th2(1:N2))
   cphi(1:N3)=cos(phi(1:N3)); sphi(1:N3)=sin(phi(1:N3))
   
   r(3)=lr1; r(6)=lr2
   do i1 = 1, N2
      do i2 = 1, N1
         do i3 = 1, N3
            call myjcb_Ne4(lr1, lr2, BR,cth1(i2),sth1(i2),cth2(i1),sth2(i1),  &
                             cphi(i3),sphi(i3),r(1),r(4),r(5),r(2))            
            call ne4pot(r, V(i3, i2, i1)) 
         end do
      end do
   end do

end 

double precision function potJA4Non(BR, lr1, lr2, theta1, theta2, phi)
   implicit none
   double precision, intent(IN) :: BR, lr1, lr2, theta1, theta2, phi

   call vNe4(lr1, lr2, BR, theta1, theta2, phi, potJA4Non)
end

double precision function pot(BR, lr1, lr2, theta1, theta2, phi)
   implicit none
   double precision, intent(IN) :: BR, lr1, lr2, theta1, theta2, phi
   double precision :: potJA4Non

   call vNe4(lr1, lr2, BR, theta1, theta2, phi, pot)
end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine vNe4(r1,r2,r0,th1,th2,phij,v)   
      implicit none
      double precision, intent(IN) :: r1, r2,r0,th1,th2, phij
      double precision, intent(OUT):: v

      double precision :: r(6)

      call myjcb_ne4(r1,r2,r0,th1,th2,phij,r(1),r(4),r(5),r(2))
      r(3) = r1 ;  r(6) = r2
      call ne4pot(r, v)
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Subroutine to calculate Ne4 potential         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine ne4pot(r, pot)
   implicit none
   double precision, intent(IN) :: r(6)
   double precision, intent(OUT):: pot

   double precision, parameter :: rmin = 0.5D0;
   double precision :: r_6;  ! r**6
   integer :: i
   
   pot = 0.0;
   do i = 1, 6
      if (r(i) < rmin)  then
         r_6 = rmin**6
      else
         r_6 = r(i)**6
      end if

      r_6 = 1.0D0/r_6
      pot = pot + r_6*r_6-r_6;

   end do

   pot = 4.0D0*pot
end 


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       b                         d
!c        \                       /
!c         \                     /
!c          e-------------------f---------------g
!c           \                 /
!c            a---------------c-------------h
!c
!c   a, b = Ne, c, d = Ne
!c   e is the COM of a and b, f is the COM of c and d:
!c
!c       the Descartes coordinates is chosen as
!c       the coordinate origin is at e(0,0,0), x axis is from e to g
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine myjcb_ne4(ab,cd,ef,c_beg,s_beg,c_dfg,s_dfg,c_befd,s_befd, ac,bc,ad,bd)      
        implicit none
        double precision, intent(IN) :: ab, cd, ef
        double precision, intent(IN) :: c_beg,s_beg, c_dfg,s_dfg,c_befd,s_befd
        double precision, intent(OUT):: ac, bc, ad, bd

        double precision :: a,b,c,d
        double precision :: ae,be,cf,df,xf,yf,zf
        double precision :: xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd
        double precision :: xac,yac,zac, xad,yad, zad
        double precision :: xbc,ybc,zbc, xbd, ybd,zbd        

        a = 20.0d0;   b = a
        c = 20.0d0;   d = c

        ae=ab*b/(a+b);     be=ab-ae; 
        cf=cd*d/(c+d);     df=cd-cf

        xf=ef;    yf=0.d0;    zf=0.d0
        xa=-ae*c_beg;   ya=-ae*s_beg;   za=0.d0
        xb= be*c_beg;   yb= be*s_beg;   zb=0.d0
        xc=ef-cf*c_dfg; yc=-cf*s_dfg*c_befd
        zc=-cf*s_dfg*s_befd

        xd=ef+df*c_dfg; yd= df*s_dfg*c_befd
        zd= df*s_dfg*s_befd

        xac=xc-xa; yac=yc-ya; zac=zc-za
        ac=dsqrt(dabs(xac*xac+yac*yac+zac*zac))  

        xad=xd-xa; yad=yd-ya; zad=zd-za
        ad=dsqrt(dabs(xad*xad+yad*yad+zad*zad))  

        xbc=xc-xb; ybc=yc-yb; zbc=zc-zb
        bc=dsqrt(dabs(xbc*xbc+ybc*ybc+zbc*zbc))  

        xbd=xd-xb; ybd=yd-yb; zbd=zd-zb
        bd=dsqrt(dabs(xbd*xbd+ybd*ybd+zbd*zbd))    
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Calculate 1D potential using fitting functions    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           1D potential for Vlr in Ne4               c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine fitVlr1(N, R1, VR1)
   implicit none
   integer, intent(IN) :: N
   double precision, intent(IN) :: R1(N)
   double precision, intent(OUT):: VR1(N)

   double precision, parameter :: coeff(12) =  (/                     &
         -48.80698836230245D0,-6.08734576622792D0, 1.99874598034173D0,& 
          52.80683314404432D0, 1.00596778725101D0, 3.89528586956230D0,&
          2.02696645335764D0, -1.85304424168344D0,-2.14057293742726D0,&
          0.06886855302831D0, -3.53913701642550D0, 2.33349809856400D0 /)

   double precision :: x0(N),x1(N)
   
   x0(1:N) = R1(1:N)**(-12)-R1(1:N)**(-6);
   call mystep(coeff(1:4), N, R1, x1)
   VR1(1:N) = x1(1:N)*x0(1:N) ;
    
   x0(1:N) = exp(0.1**(6)*R1(1:N)**(-6))+coeff(10)*exp(coeff(11)*(R1(1:N)-coeff(12))**(2));
   call mystep(coeff(5:8), N, R1, x1)
    
   VR1(1:N) = VR1(1:N) + x0(1:N)*x1(1:N) + coeff(9);

end 

!***************************************************************
subroutine fitVlr2(N, R2, VR2)
   implicit none
   integer, intent(IN) :: N
   double precision, intent(IN) :: R2(N)
   double precision, intent(OUT):: VR2(N)

   call fitVlr1(N, R2, VR2)
end 

!************************************************************
!************************************************************
!************************************************************nd

!************************************************************
subroutine fitVBR(N, R3, VR3)
   implicit none
   integer, intent(IN) :: N
   double precision, intent(IN) :: R3(N)
   double precision, intent(OUT):: VR3(N)

   double precision, parameter :: coeff(33) = (/                    &
      -31.14406978428782D0, 7.08391290586332D0, 2.20301370051600D0, &
      -11.22846955152019D0, 3.40601563427331D0, 8.20761388563837D0, &
        0.78516268924461D0,-3.68814729848298D0, 1.38449071602294D0, &
       -7.82646131752593D0, 1.61128697660778D0, 6.46576340978716D0, &
        1.06413961335373D0,-2.28331753072536D0, 2.57762485755913D0, &
       16.76245458588856D0, 8.81711575699348D0,-7.58229455933339D0, &
        1.71655615718333D0, 3.14257087001422D0,16.15662966048344D0, &
       -0.53061441313732D0,-0.39260653137463D0,-0.99417317403920D0, &
        2.23164995612228D0, 0.91793288482833D0, 0.78542452000807D0, &
        2.45878772643765D0, 2.22726024780452D0, 1.86336101108059D0, &
        0.39705056650805D0, 0.33108157531040D0, 0.06077136846554D0 /)
   
   double precision :: x1(N),x2(N),x(3)

   x1(1:N) = exp(coeff(22)*R3(1:N)**(2)+coeff(23)*R3(1:N)+coeff(24)) 
   call mystep(coeff(1:4), N, R3, x2)
   VR3(1:N) = x2(1:N)*x1(1:N)

   x1 = coeff(25)*(1.0D0-exp(-coeff(31)*(R3(1:N)-coeff(26)))**(2))
   call mymidFilter(coeff(5:10), N, R3, x2)
   VR3(1:N) = VR3(1:N) + x2(1:N)*x1(1:N)

   x1 = coeff(27)*(1.0D0-exp(-coeff(32)*(R3(1:N)-coeff(28)))**(2))   
   call mymidFilter(coeff(11:16), N, R3, x2)
   VR3(1:N) = VR3(1:N) + x2(1:N)*x1(1:N)
   
   x1 = coeff(29)*(1.0D0-exp(-coeff(33)*(R3(1:N)-coeff(30)))**(2))
   call mystep(coeff(17:20), N, R3, x2)
   VR3(1:N) = VR3(1:N) + x2(1:N)*x1(1:N) + coeff(21)  
  
end 


