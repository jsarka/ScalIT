!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!c  Compex version of QMR part of OSB:                              c
!c  Subroutines:                                                    c
!c    integer function OSB_QMRDX(N,B,X, ErRR)                       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      General Case of  X=B/(H-E-iAP)                              c
!c            P=1/(EigVal-E-iAP)                                    c
!c  Normally: sHC=sPC.  The user needs to call initOSBW(0/1) before c
!c   call this function and call finalOSBW() after finish work for  c
!c              OSBW preconditioner                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function OSB_QMRDX(NIN, B, X, ERR_RES)
      integer, intent(IN)   :: NIN
      double complex, intent(in)  :: B(NIN)
      double complex, intent(out) :: X(NIN)
      double precision, intent(out)   :: ERR_RES
 
      integer :: qmr_cx
 
      OSB_QMRDX = QMR_cx(sQMRConvType,sQMR%mMax,sQMR%mTol,NIN,   &
                    B, HX_DX, PX_DX, X, ERR_RES)     

  end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   X=B/(H+iAP-E), P=1/(EigVal-E) or P=1/(EigVal+iAP-E)      c
!c        Preconditioner is determined by sPC                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMRDX_EA(NIN, B, X, ERR_RES)
      integer, intent(IN)   :: NIN
      double complex, intent(in)  :: B(NIN)
      double complex, intent(out) :: X(NIN)
      double precision, intent(out)   :: ERR_RES

      integer :: qmr_cx

      QMRDX_EA = QMR_cx(sQMRConvType, sQMR%mMax,sQMR%mTol,NIN,  &
                    B, EHX_DX, EPX_DX, X, ERR_RES) 

  end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMRDX_PA(NIN, B, X, ERR_RES)
      integer, intent(IN)   :: NIN
      double complex, intent(in)  :: B(NIN)
      double complex, intent(out) :: X(NIN)
      double precision, intent(out)   :: ERR_RES

      integer :: qmr_cx

      QMRDX_PA = QMR_cx(sQMRConvType, sQMR%mMax,sQMR%mTol,NIN,  &
                    B, EHX_PDX, EPX_PCX, X, ERR_RES) 

  end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   X=B/(H-iAP-E), P=1/(EigVal-E) or P=1/(EigVal-iAP-E)      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMRDX_MA(NIN, B, X, ERR_RES)
      integer, intent(IN)   :: NIN
      double complex, intent(in)  :: B(NIN)
      double complex, intent(out) :: X(NIN)
      double precision, intent(out)   :: ERR_RES

      integer :: qmr_cx

      QMRDX_MA = QMR_cx(sQMRConvType, sQMR%mMax, sQMR%mTol, NIN, &
                    B, EHX_MDX, EPX_MCX, X, ERR_RES) 

  end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

