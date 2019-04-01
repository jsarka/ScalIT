!
! Module to get wavefunction for tri-atomic molecules
! using Jacobi coorinates
!

module wave3jb
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

end module


