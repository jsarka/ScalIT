!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Initialize and update  VOUT                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine VINIT(NOUT, VOUT)
      implicit none

      integer,intent(IN) ::  NOUT                
      double precision, intent(OUT):: VOUT(NOUT, NOUT)  

      integer   ::  i

      Vout(1:NOUT, 1:Nout) = 0.0D0

      do i = 1,NOUT
         VOUT(i,i) = 1.0D0
      end do

  end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine VUPDATE(P, Q, PHI, NOUT, VOUT)
      implicit none
      integer,intent(IN) ::  P, Q ,NOUT
      double precision, intent(IN)   :: PHI 
      double precision, intent(INOUT)::VOUT(NOUT, NOUT) 

      double precision :: tmp(NOUT) 
      double precision :: C, S  

      C = DCOS(PHI);   S = DSIN(PHI)

      tmp(1:NOUT)     = C * VOUT(1:NOUT,P) - S * VOUT(1:NOUT, Q)
      VOUT(1:NOUT, Q) = C * VOUT(1:NOUT,Q) + S * VOUT(1:NOUT, P)
      VOUT(1:NOUT, P) = tmp(1:NOUT)

  end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


