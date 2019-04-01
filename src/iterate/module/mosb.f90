!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Module dealing with communication directly      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
module mosb
     use mosbtype
     implicit none
!     private
     include 'mosb.data.h'
     include 'mosbw.data.h'  
     include 'mosb3.data.h'

  contains
     include 'mosb.util.h'
     include 'mosb.init.h'
     include 'mosb3.init.h'

     include 'mosb.io.h'
     include 'mosb.print.h'
!     include 'mosb.io1.h'
     include 'mosb.io10.h'

     include 'mosb.diaginit.h'
     include 'mosb.diag.h'
     include 'mosb.progDiag.h'
     include 'mosb3.updatex.h'
     
     include 'mosb.hx.h'
     include 'mosb3.hx.h'
     include 'mosb3.mxgrid.h'
     include 'mosb3.mxseq.h'

     include 'mosb.hxcx.h'
     include 'mosb3.hxcx.h'
     include 'mosb3.mxgrid_cx.h'
     include 'mosb3.mxseq_cx.h'

     include 'mosb.hxdx.h'
     include 'mosb3.hxdx.h'
     include 'mosb3.mxgrid_dx.h'
     include 'mosb3.mxseq_dx.h'

     include 'mosb.progHX.h'

     include 'mosb.px.h'
     include 'mosb3.px.h'
     include 'mosb3.vxgrid.h'
     include 'mosb3.vxseq.h'

     include 'mosb.pxdx.h'
     include 'mosb3.pxdx.h'
     include 'mosb3.vxgrid_dx.h'
     include 'mosb3.vxseq_dx.h'

     include 'mosb.progPX.h'

!cccccccccc OSBW cccccccccccccccccc
     include 'mosbw.init.h'
     include 'mosbw.io.h'
     include 'mosbw.osbw.h'
     include 'mosbw.hmat.h'
     include 'mosbw.hij.h'
     include 'mosbw.hij1.h'
     include 'mosbw.hijcx.h'
     include 'mosbw.hijcx1.h'
     include 'mosbw.vi.h'
     include 'mosbw.vx.h'
     include 'mosbw.pij.h'
     include 'mosbw.progHij.h'
     !include 'mosbw.progPij.h'

end module
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
