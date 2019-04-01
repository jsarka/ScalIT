! calculate wave functions with local modes:
! 
!
module wave3lm

   implicit none
   private

   interface readW3
      module procedure readW3StdIO
      module procedure readW3File
   end interface
   
   public :: initW3, finalW3, printW3, calSaveWF

   include "wave3.data"

 contains

   include "wave3.init"
   include "wave3.io"
   include "wave3jb.init"
   include "wave3jb.io"
   include "wave3jb.wfp1"
   include "wave3jb.wfm1"
   include "wave3jb.wfs1"
   include "wave3jb.wfp2"
   include "wave3jb.wfm2"
   include "wave3jb.wfs2"

!
! Switch internal coordinates: 
! (A - [B) - C] : internal coordinates(r1, r2, cos(a)) 
! A-(BC) : Jacobi coordinates, lr, BR, cos(theta))
!
   subroutine switchJB(r1, r2, ca, lr, BR, cth)
      double precision, intent(IN)  :: r1, r2, ca
      double precision, intent(OUT) :: lr, BR, cth
    
      double precision :: rB

      lr = r2; rB = mass(3)/(mass(2)+mass(3)); 
      BR = r1**2+rB**2-2*rB*r1*ca
      BR = sqrt(BR)
      if ((rB==0.0D0) .OR. (BR==0.0D0)) then
         cth=1.0D0
      else
         cth = (rB**2+BR**2-r1**2)/(2.0*rB*BR)
         if (cth>1.0D0) then
            cth=1.0D0
         else
            if (cth<-1.0D0) cth=-1.0D0
         end if
      end if
   end subroutine
end module

