!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Module dealing with communication directly      c
!c     This is MPI version without parallel IO         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
module posb
     use mosbtype
     implicit none
     include 'mosb.data.h'
     include 'mosbw.data.h'  
     include 'mosb3.data.h'

  contains
     include 'mosb.util.h'
     include 'mosb.init.h'
     include 'mosb3.init.h'

     include 'posb.io.h'
     include 'posb.io1.h'
     include 'posb.io2.h'
     include 'mosb.print.h'


     include 'mosb.diaginit.h'
     include 'posb.diag.h'
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

     include 'posb.progHX.h'

     include 'mosb.px.h'
     include 'mosb3.px.h'
     include 'mosb3.vxgrid.h'
     include 'mosb3.vxseq.h'

     include 'mosb.pxdx.h'
     include 'mosb3.pxdx.h'
     include 'mosb3.vxgrid_dx.h'
     include 'mosb3.vxseq_dx.h'

     include 'posb.progPX.h'

!cccccccccc OSBW cccccccccccccccccc
     include 'mosbw.init.h'
     include 'mosbw.io.h'
     include 'mosbw.osbw.h'
     include 'mosbw.hmat.h'
     include 'posbw.hij.h'
     include 'posbw.hij1.h'
     include 'posbw.hijcx.h'
     include 'posbw.hijcx1.h'
     include 'mosbw.vi.h'
     include 'mosbw.vx.h'
     include 'mosbw.pij.h'
     include 'mosbw.progHij.h'

end module
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
