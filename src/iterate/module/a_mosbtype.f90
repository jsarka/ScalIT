!cccccccccccccccccccccccccccccccccccccccccc
!c     Data structures used in MOSB       c
!cccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccc
module mosbType

   implicit none
   include 'mpif.h'
   integer, parameter :: FMAX = 60   
   include 'mosbtype.data.h'

   contains

     include 'mosbtype.print.h'
     include 'mosbtype.cal.h'
     include 'mosbtype.mosb.h'

end module
!cccccccccccccccccccccccccccccccccccccccccc

