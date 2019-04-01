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
  integer function OSB_QMRDX(N, B, X, ERR_RES)
      integer, intent(IN)   :: N
      double complex, intent(in)  :: B(N)
      double complex, intent(out) :: X(N)
      double precision, intent(out)   :: ERR_RES
 
      integer :: qmr_cx_mpi
 
      OSB_QMRDX = QMR_CX_MPI(sQMRConvType,sQMR%mMax,sQMR%mTol,N,   &
                    B, HX_DX, PX_DX, X, ERR_RES)     

  end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   X=B/(H+iAP-E), P=1/(EigVal-E) or P=1/(EigVal+iAP-E)      c
!c        Preconditioner is determined by sPC                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMRDX_EA(N, B, X, ERR_RES)
      integer, intent(IN)   :: N
      double complex, intent(in)  :: B(N)
      double complex, intent(out) :: X(N)
      double precision, intent(out)   :: ERR_RES

      integer :: qmr_cx_mpi

      QMRDX_EA = QMR_CX_MPI(sQMRConvType, sQMR%mMax,sQMR%mTol,N,  &
                    B, EHX_DX, EPX_DX, X, ERR_RES) 

  end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   X=B/(H-iAP-E), P=1/(EigVal-E) or P=1/(EigVal-iAP-E)      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMRDX_PA(N, B, X, ERR_RES)
      integer, intent(IN)   :: N
      double complex, intent(in)  :: B(N)
      double complex, intent(out) :: X(N)
      double precision, intent(out)   :: ERR_RES

      integer :: qmr_cx_mpi

      QMRDX_PA = QMR_CX_MPI(sQMRConvType, sQMR%mMax, sQMR%mTol, N, &
                    B, EHX_PDX, EPX_PCX, X, ERR_RES) 

  end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   X=B/(H-iAP-E), P=1/(EigVal-E) or P=1/(EigVal-iAP-E)      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMRDX_MA(N, B, X, ERR_RES)
      integer, intent(IN)   :: N
      double complex, intent(in)  :: B(N)
      double complex, intent(out) :: X(N)
      double precision, intent(out)   :: ERR_RES

      integer :: qmr_cx_mpi

      QMRDX_MA = QMR_CX_MPI(sQMRConvType, sQMR%mMax, sQMR%mTol, N, &
                    B, EHX_MDX, EPX_MCX, X, ERR_RES) 

  end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
