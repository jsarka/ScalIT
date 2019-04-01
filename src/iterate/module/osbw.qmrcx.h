!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!c  Compex version of QMR part of OSB:                              c
!c  Subroutines:                                                    c
!c  integer function OSB_QMRCX_PX(N,B,X, ErRR)                      c
!c  integer function OSB_QMRCX_MX(N,B,X, ErRR)                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      General Case of  X=B/(H-E[+/-iAP,iAPP,iAPR]),               c
!c            P=1/(EigVal-E [+/-iAP,iAPP,iAPR])                     c
!c  Normally: sHC=sPC.  The user needs to call initOSBW(0/1) before c
!c   call this function and call finalOSBW() after finish work for  c
!c              OSBW preconditioner                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function OSB_QMRCX(NIN, B, X, ERR_RES)
      integer, intent(IN)   :: NIN
      double complex, intent(in)  :: B(NIN)
      double complex, intent(out) :: X(NIN)
      double precision, intent(out)   :: ERR_RES
 
      integer :: qmr_cx
 
      OSB_QMRCX = QMR_cx(sQMRConvType,sQMR%mMax,sQMR%mTol,NIN,   &
                    B, HX_CX, PX_DX, X, ERR_RES)     

  end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   X=B/(H+iAP-E), P=1/(EigVal-E) or P=1/(EigVal+iAP-E)      c
!c        Preconditioner is determined by sPC                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMRCX_EA(NIN, B, X, ERR_RES)
      integer, intent(IN)   :: NIN
      double complex, intent(in)  :: B(NIN)
      double complex, intent(out) :: X(NIN)
      double precision, intent(out)   :: ERR_RES

      integer :: qmr_cx

      QMRCX_EA = QMR_cx(sQMRConvType, sQMR%mMax,sQMR%mTol,NIN,  &
                    B, EHX_PCX, EPX_DX, X, ERR_RES) 

  end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   X=B/(H+iAP-E), P=1/(EigVal-E) or P=1/(EigVal+iAP-E)      c
!c        Preconditioner is determined by sPC                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMRCX_PA(NIN, B, X, ERR_RES)
      integer, intent(IN)   :: NIN
      double complex, intent(in)  :: B(NIN)
      double complex, intent(out) :: X(NIN)
      double precision, intent(out)   :: ERR_RES

      integer :: qmr_cx

      QMRCX_PA = QMR_cx(sQMRConvType, sQMR%mMax,sQMR%mTol,NIN,  &
                    B, EHX_PCX, EPX_PCX, X, ERR_RES) 

  end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   X=B/(H-iAP-E), P=1/(EigVal-E) or P=1/(EigVal-iAP-E)      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMRCX_MA(NIN, B, X, ERR_RES)
      integer, intent(IN)   :: NIN
      double complex, intent(in)  :: B(NIN)
      double complex, intent(out) :: X(NIN)
      double precision, intent(out)   :: ERR_RES

      integer :: qmr_cx

      QMRCX_MA = QMR_cx(sQMRConvType, sQMR%mMax, sQMR%mTol, NIN, &
                    B, EHX_MCX, EPX_MCX, X, ERR_RES) 

  end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

