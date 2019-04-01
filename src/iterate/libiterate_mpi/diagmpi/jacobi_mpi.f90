!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C                           MJACOBI function                                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  double precision function MJACOBI_MPI(COMM,P,Q,NOUT,NIN,HOUT,IERR)
      implicit none
      include 'mpif.h'
      integer,  intent(IN)  :: COMM, P, Q, NOUT,NIN               
      double precision,intent(INOUT)::HOUT(NIN,NOUT,NOUT) 
      integer,  intent(OUT) :: IERR

      double precision, parameter :: PIE = 3.14159265358979323846264338328d0

      integer ::  R                       
      double precision :: A, B, C, X       
      double precision :: SN, CS, C2, S2, SC  
        
      double precision :: DDG(NIN), ODG(NIN)
      double precision :: NEWP(NIN), NEWQ(NIN), NEWQP(NIN)
      double precision :: LOCABC(3), GLBABC(3)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Calculate rotation angle. It's weird that using XSUBYVEC 
!c     replace AXBYVEC gives the worse results, so use AXBYVEC
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ODG(1:NIN) = HOUT(1:NIN, P, Q)
      DDG(1:NIN) = 0.5D0*HOUT(1:NIN, P, P) - 0.5D0*HOUT(1:NIN, Q, Q)

      LOCABC(1) = dot_product(ODG(1:NIN),ODG(1:NIN))
      LOCABC(2) = 2.0D0 * dot_product(ODG(1:NIN),DDG(1:NIN))
      LOCABC(3) = dot_product(DDG(1:NIN),DDG(1:NIN))

      call MPI_ALLREDUCE(LOCABC, GLBABC, 3, MPI_DOUBLE_PRECISION,  &
                  MPI_SUM, COMM,IERR )
      A = GLBABC(1) ;      B = GLBABC(2) ;      C = GLBABC(3) 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Determine rotation angle. NOTE: Exception handling a bit different
!c     than before; fixes a slight inconsistency in the previous version.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (A.eq.C) then
         if (B.eq.0d0) then
            X = 0.0d0
         else 
            X = -DSIGN(PIE/8.0d0,B)
         end if 
      else
         X = DATAN(B/(A-C))/4.0d0
         if (A.gt.C) then
            if (X.lt.0) then
               X = X + PIE/4.0d0
            else
               X = X - PIE/4.0d0
            end if
         end if
      end if

      SN = DSIN(X);    CS = DCOS(X);
      C2 = CS**2;      S2 = SN**2;      SC = SN*CS

!cccccccccccccccccccccccccccccccccccccccc

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
          if (R.ne.P.and.R.ne.Q) then 
             NEWP(1:NIN) = CS*HOUT(1:NIN,R,P) - SN*HOUT(1:NIN,R,Q)             
             NEWQ(1:NIN) = CS*HOUT(1:NIN,R,Q) + SN*HOUT(1:NIN,R,P)             
             HOUT(1:NIN,R,P) = NEWP(1:NIN)
             HOUT(1:NIN,P,R) = NEWP(1:NIN)
             HOUT(1:NIN,R,Q) = NEWQ(1:NIN)
             HOUT(1:NIN,Q,R) = NEWQ(1:NIN)
         end if
      end do
         
      MJACOBI_MPI = X 

  end

