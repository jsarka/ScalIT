!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 double precision   function  MJACOBI(P,Q,NOUT,NIN,HOUT)
      implicit none
      integer,intent(IN) :: P, Q, NOUT,NIN              
      double precision, intent(INOUT)::  HOUT(NIN, NOUT, NOUT)   

      double precision,parameter :: PIE = 3.14159265358979323846264338328d0

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision  ::  A, B, C, X  
      double precision  ::  SN, CS, C2, S2, SC 

      integer  :: R                           

      double precision :: DDG(NIN), ODG(NIN) 
      double precision :: NEWP(NIN), NEWQ(NIN), NEWQP(NIN)   

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Calculate rotation angle. It's weird that using XSUBYVEC  c
!c     replace AXBYVEC gives the worse results, so use AXBYVEC   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ODG(1:NIN) = HOUT(1:NIN, P, Q)
      DDG(1:NIN) = 0.5D0*HOUT(1:NIN, P, P) - 0.5D0*HOUT(1:NIN, Q, Q)

      A = dot_product(ODG(1:NIN), ODG(1:NIN))        
      B = 2.0D0 * dot_product(ODG(1:NIN), DDG(1:NIN)) 
      C = dot_product(DDG(1:NIN), DDG(1:NIN))         

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Determine rotation angle. NOTE: Exception handling a bit different   c
!c   than before; fixes a slight inconsistency in the previous version.   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (A == C) then
         if (B == 0d0) then
            X = 0.0d0
         else 
            X = -DSIGN(PIE/8.0d0,B)
         end if 
      else
         X = DATAN(B/(A-C))/4.0d0
         if (A > C) then
            if (X  < 0.0D0 ) then
               X = X + PIE/4.0d0
            else
               X = X - PIE/4.0d0
            end if
         end if
      end if


      SN = DSIN(X);  CS = DCOS(X)
      C2 = CS**2;    S2 = SN**2;      SC = SN*CS

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Update HOUT matrix, set intersection blocks first      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      NEWP(1:NIN) = C2*HOUT(1:NIN,P,P) + S2*HOUT(1:NIN,Q,Q) - 2*SC*ODG(1:NIN)
      NEWQ(1:NIN) = C2*HOUT(1:NIN,Q,Q) + S2*HOUT(1:NIN,P,P) + 2*SC*ODG(1:NIN)
      NEWQP(1:NIN)= (C2-S2)*ODG(1:NIN) + 2*SC*DDG(1:NIN)
      HOUT(1:NIN,P,P) = NEWP(1:NIN)
      HOUT(1:NIN,Q,Q) = NEWQ(1:NIN)
      HOUT(1:NIN,P,Q) = NEWQP(1:NIN)
      HOUT(1:NIN,Q,P) = NEWQP(1:NIN)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
!c     set p, q row/column blocks; sign change from before   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do R=1,NOUT
         if ( (R /= P).AND. (R /= Q)) then 
            NEWP(1:NIN) = CS*HOUT(1:NIN,R,P) - SN*HOUT(1:NIN,R,Q)             
            NEWQ(1:NIN) = CS*HOUT(1:NIN,R,Q) + SN*HOUT(1:NIN,R,P)             
            HOUT(1:NIN,R,P) = NEWP(1:NIN);  HOUT(1:NIN,P,R) = NEWP(1:NIN)
            HOUT(1:NIN,R,Q) = NEWQ(1:NIN);   HOUT(1:NIN,Q,R) = NEWQ(1:NIN)
         end if
      end do
         
      MJACOBI = X      

end
