!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!c  QMR part of OSB:                                              c
!c  Subroutines:                                                  c
!c  QMR_OSB0 (N, B, X, ERES)  ! General (sHC,sPC) preconditioner  c
!c  QMR_OSB  (N, B, X, ERES)  ! OSB preconditioner                c
!c  QMR_OSBD1(N, B, X, ERES)  ! OSBD1 preconditioner              c
!c  QMR_OSBD2(N, B, X, ERRS)  ! OSBD2 preconditioner              c
!c  QMR_OSBW (N, B, X, ERES)  ! OSBW preconditioner               c
!c        ! use HOSB if it is stored in memory(10 times faster)   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function OSB_QMR(NIN, B, X, ERES)
      integer, intent(IN) :: NIN
      double precision, intent(in) :: B(NIN) 
      double precision, intent(out) :: X(NIN), ERES

      integer :: qmr

      OSB_QMR = QMR(sQMRConvType, sQMR%mMax, sQMR%mTol, NIN,     &
                    B, HX, PX, X, ERES)      

  end function
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMR_OSB(NIN, B, X, ERES)
      integer, intent(IN) :: NIN
      double precision, intent(in) :: B(NIN) 
      double precision, intent(out) :: X(NIN), ERES

      integer :: qmr

      QMR_OSB = QMR(sQMRConvType, sQMR%mMax, sQMR%mTol, NIN,     &
                    B, EHX, EPX, X, ERES)

  end function
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMR_OSB0(NIN, B, X, ERES)
      integer, intent(IN) :: NIN
      double precision, intent(in) :: B(NIN) 
      double precision, intent(out) :: X(NIN), ERES

      integer :: qmr

      QMR_OSB0 = QMR(sQMRConvType, sQMR%mMax, sQMR%mTol, NIN,    &
                    B, EHX, EPX0, X, ERES)

  end function
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMR_OSBD1(NIN, B, X, ERES)
      integer, intent(IN) :: NIN
      double precision, intent(in) :: B(NIN) 
      double precision, intent(out) :: X(NIN), ERES

      integer :: qmr

      QMR_OSBD1 = QMR(sQMRConvType, sQMR%mMax, sQMR%mTol, NIN,   &
                     B, EHX, EPXD1, X, ERES)

  end function 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMR_OSBD2(NIN, B, X, ERES)
      integer, intent(IN) :: NIN
      double precision, intent(in) :: B(NIN) 
      double precision, intent(out) :: X(NIN), ERES

      integer :: qmr

      QMR_OSBD2 = QMR(sQMRCONVTYPE, sQMR%mMAX, sQMR%mTol, NIN,   &
                     B, EHX, EPXD2, X, ERES)

  end function 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c          Perform QMR using the first order of Hij              c
!c  All first order Hij within the energy window are calculated   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    initOSBW should be called before calling this function      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMR_OSBW(NIN, B,  X, ERES)  
      integer, intent(IN) :: NIN
      double precision, intent(in) :: B(NIN) 
      double precision, intent(out) :: X(NIN), ERES
      
      integer :: qmr
  
      QMR_OSBW = QMR(sQMRCONVTYPE,sQMR%mMAX,sQMR%mTOL,NIN,       &
               B, EHX, EPXW, X, ERES)

  end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c initOSBW not need to be called before calling this function    c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMR_OSBW_NOMEM(NIN, B,  X, ERES)  
      integer, intent(IN) :: NIN
      double precision, intent(in) :: B(NIN)
      double precision, intent(out) :: X(NIN), ERES
      
      integer :: qmr
  
      if ( INITOSBW()) then
          QMR_OSBW_NOMEM=QMR(sQMRCONVTYPE,sQMR%mMAX,sQMR%mTOL,NIN,&
                         B, EHX, EPXW, X, ERES)
      else
          QMR_OSBW_NOMEM = 0
      end if
 
      call FINALOSBW()

  end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
