!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    Calculate JB parameters when exchange atoms in A3        c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine jbech(Br, lr, cth, BR1, lr1, cth1, BR2, lr2, cth2)
   implicit none
   double precision, intent(IN)  :: BR, lr, cth
   double precision, intent(OUT) :: BR1,lr1,cth1,BR2,lr2,cth2

   double precision :: r12, r22

   lr1 = BR**2 + (0.5D0*lr)**2 - BR*lr*cth
   lr2 = BR**2 + (0.5D0*lr)**2 + BR*lr*cth

   BR1 = DSQRT(0.5D0*(lr**2+lr2-0.5D0*lr1))
   BR2 = DSQRT(0.5D0*(lr**2+lr1-0.5D0*lr2))

   lr1 = DSQRT(lr1);   lr2 = DSQRT(lr2)
   cth1 = (BR1**2 + (lr1/2.0D0)**2 - lr2**2)/(BR1*lr1)
   cth2 = (BR2**2 + (lr2/2.0D0)**2 - lr**2)/(BR2*lr2)

   if (cth1<-1.0D0) then
      cth1=-1.0D0
   else
      if (cth1>1.0D0) cth1=1.0D0
   end if

   if (cth2<-1.0D0) then
      cth2=-1.0D0
   else
      if (cth2>1.0D0) cth2=1.0D0
   end if

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Transformation between Jacobi and Hyperspherical cooridnates    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Jacobi coorinate( R0/R, r1/r, phi(angle)) =>                   c
!c  Hyperspherical coordinate(pho,tantheta, cos2chi, sin2chi)      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine jb2hs(BS, ss, phi, rho, theta, chi )
   implicit none
   double precision, intent(IN) :: BS, ss, phi
   double precision, intent(OUT):: rho, theta, chi

   double precision :: rd2, sr0, cr0, rm 

   rho = sqrt(BS**2 + ss**2)

   rd2 = BS**2 - ss**2
   sr0 = 2.0D0*BS*ss*sin(phi)
   cr0 = 2.0D0*BS*ss*cos(phi)
   rm = sqrt(rd2**2+cr0**2)

   if (sr0==0.0D0) then
      theta = 0.5D0*ACOS(-1.0D0)
   else
      theta = ATAN( rm / sr0 )
   end if

   if (rm==0.0D0) then
      chi = 0.0D0
   else
      chi = rd2/rm
      if (chi>1.0D0) chi=1.0D0
      if (chi<-1.0D0) chi=-1.0D0
      chi = ACOS(chi)
      if (sr0 < 0.0D0) chi = 2.0D0*acos(-1.0D0)-chi
   end if 
   chi = 0.5D0*chi
end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Hyperspherical coorinate( pho, theta, chi) =>                  c
!c      Jacobi coordinate(R0, r1, cphi)                            c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine hs2jb(rho, theta, chi, BS, ss, cphi)
   implicit none
   double precision, intent(IN) :: rho, theta, chi
   double precision, intent(OUT):: BS, ss, cphi

   double precision, parameter :: SQRT1_2=0.7071067811865475244D0  !1/sqrt(2)

   double precision :: a0, a1, a2, s2chi, c2chi, sth

   a0=rho*SQRT1_2;        sth=sin(theta)
   s2chi=sin(2.0D0*chi);  c2chi=cos(2.0D0*chi)
   a1=sth*c2chi;          a2=sth*s2chi     
   BS=a0*DSQRT(1.0D0+a1)
   ss=a0*DSQRT(1.0D0-a1)

   a1=a1**2
   if (a1==1.0D0) then
       if (a2>0.0D0) then
          cphi = 1.0D0
       else
          cphi = -1.0D0
       end if
   else
       cphi = sth*s2chi/DSQRT(1.0D0-a1)
   end if

   if (cphi<-1.0D0) then
      cphi=-1.0D0
   else
      if (cphi>1.0D0) cphi=1.0D0
   end if
end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine xy2jb(rho, x0, y0, BS, ls, cth, F0)
   implicit none
   double precision, intent(IN)  :: rho, x0, y0
   double precision, intent(OUT) :: BS, ls, cth, F0
   double precision, parameter :: SQRT1_2=0.7071067811865475244D0  !1/sqrt(2)
   double precision :: xy2, stheta, ctheta, s2chi, c2chi

   xy2 = x0**2 + y0**2
   if (xy2==0.0D0) then
      BS=0.0D0; ls=0.0D0; cth=0.0D0; F0=0.0D0
   else
      stheta=2.0D0*DSQRT(xy2)/(1.0D0+xy2)
      ctheta=(1.0D0-xy2)/(1+xy2)
      s2chi=2*x0*y0/xy2
      c2chi=(x0**2-y0**2)/xy2
      xy2=stheta*c2chi
      BS=SQRT1_2*rho*DSQRT(1.0D0+xy2)
      ls=SQRT1_2*rho*DSQRT(1.0D0-xy2)
      xy2=xy2**2
      if (xy2==1.0D0) then
         cth=0.0D0
      else
         cth=stheta*s2chi/DSQRT(1.0-xy2)
      end if
      F0=2.0D0*stheta*ctheta*rho**5
   end if

   if (cth<-1.0D0) then
      cth=-1.0D0
   else
      if (cth>1.0D0) cth=1.0D0
   end if
end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
