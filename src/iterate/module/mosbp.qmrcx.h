!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!c  Compex version of QMR part of OSB:                                      c
!c  Subroutines:                                                            c
!c  integer function OSB_QMRCX_PX(N,B,X, ErRR)                              c
!c  integer function OSB_QMRCX_MX(N,B,X, ErRR)                              c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      General Case of  X=B/(H-E[+/-iAP,iAPP,iAPR]),               c
!c            P=1/(EigVal-E [+/-iAP,iAPP,iAPR])                     c
!c  Normally: sHC=sPC.  The user needs to call initOSBW(0/1) before c
!c   call this function and call finalOSBW() after finish work for  c
!c              OSBW preconditioner                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function OSB_QMRCX(N, B, X, ERR_RES)
      integer, intent(IN)   :: N
      double complex, intent(in)  :: B(N)
      double complex, intent(out) :: X(N)
      double precision, intent(out)   :: ERR_RES
 
      integer :: qmr_cx_mpi
 
      OSB_QMRCX = QMR_CX_MPI(sQMRConvType,sQMR%mMax,sQMR%mTol,N,   &
                    B, HX_CX, PX_DX, X, ERR_RES)     

  end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   X=B/(H+iAP-E), P=1/(EigVal-E) or P=1/(EigVal+iAP-E)      c
!c        Preconditioner is determined by sPC                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMRCX_EA(N, B, X, ERR_RES)
      integer, intent(IN)   :: N
      double complex, intent(in)  :: B(N)
      double complex, intent(out) :: X(N)
      double precision, intent(out)   :: ERR_RES

      integer :: qmr_cx_mpi

      QMRCX_EA = QMR_CX_MPI(sQMRConvType, sQMR%mMax,sQMR%mTol,N,  &
                    B, EHX_CX, EPX_DX, X, ERR_RES) 

  end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   X=B/(H-iAP-E), P=1/(EigVal-E) or P=1/(EigVal-iAP-E)      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMRCX_PA(N, B, X, ERR_RES)
      integer, intent(IN)   :: N
      double complex, intent(in)  :: B(N)
      double complex, intent(out) :: X(N)
      double precision, intent(out)   :: ERR_RES

      integer :: qmr_cx_mpi

      QMRCX_PA = QMR_CX_MPI(sQMRConvType, sQMR%mMax, sQMR%mTol, N, &
                    B, EHX_PCX, EPX_PCX, X, ERR_RES) 

  end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   X=B/(H-iAP-E), P=1/(EigVal-E) or P=1/(EigVal-iAP-E)      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMRCX_MA(N, B, X, ERR_RES)
      integer, intent(IN)   :: N
      double complex, intent(in)  :: B(N)
      double complex, intent(out) :: X(N)
      double precision, intent(out)   :: ERR_RES

      integer :: qmr_cx_mpi

      QMRCX_MA = QMR_CX_MPI(sQMRConvType, sQMR%mMax, sQMR%mTol, N, &
                    B, EHX_MCX, EPX_MCX, X, ERR_RES) 

  end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
