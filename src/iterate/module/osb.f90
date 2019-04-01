!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Module for OSB package: Just provide basic operation   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
module OSB
   use OSBType
   implicit none

!   private
   include "osb.data.h"
   include "osb.interface.h"
   
   contains

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     logical function initOSB()
         call calOSB()
         initOSB = allocOSB()
     end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     subroutine finalOSB()
         call deallocOSB()
     end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include "osb.init.h"
      include "osb.io.h"
      include "osb.print.h"

      include "osb.diag.h"
      include "osb.diag.init.h"
      include "osb.progDiag.h"   ! Program for Block_Jacobi diag
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      include "osb.hx.h"     
      include "osb.hxdx.h"
      include "osb.hxcx.h"
      include "osb.progHX.h"     ! Testing program for H*X
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include "osb.osbw.h"
      include "osb.px.h"
      include "osb.pxdx.h"
      include "osb.mypxdx.h"
      include "osb.progPX.h"     ! Testing program for P*X
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include "osb.hij.h"
      include "osb.hijcx.h"
      include "osb.pij.h"
      include "osb.vi.h"      
      include "osb.vx.h"  
!      include "osb.vmat.h"

      ! mainly used for testing of Hij
      include "osb.vi1.h"
      include "osb.index.h"
      include "osb.hij1.h"
      include "osb.progHij.h"    ! Program to calculate Hij
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

end module
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
