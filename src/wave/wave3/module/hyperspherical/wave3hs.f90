!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Module to get wavefunction for tri-atomic molecules      c
!c            in hyperspherical coordinates                   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

module wave3hs
   implicit none
   private

   interface readW3
      module procedure readW3StdIO
      module procedure readW3File
   end interface

   public :: initW3, finalW3
   public :: printW3,calSaveWF

   include "wave3.data"

 contains
   include "wave3.init"
   include "wave3.io"
   include "wave3hs.init"
   include "wave3hs.io"
   include "wave3hs.wfp1"
   include "wave3hs.wfm1" 
   include "wave3hs.wfs1" 
   include "wave3hs.wfp2" 
   include "wave3hs.wfm2" 
   include "wave3hs.wfs2" 

end module

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

