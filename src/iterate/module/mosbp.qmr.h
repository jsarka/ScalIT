!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!c  QMR: linear solver using OSB, OSBD, OSBW preconditioner       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function OSB_QMR( N, B, X, ERR_RES)
      integer, intent(IN)     :: N
      double precision, intent(in)  :: B(N) 
      double precision, intent(out) :: X(N) 
      double precision, intent(out) :: ERR_RES

      integer :: qmr_MPI

      !if (id==rootID) print *,'Iter:', sQMR%mMax,sQMR%mTol

      OSB_QMR = QMR_MPI(sQMRConvType, sQMR%mMax, sQMR%mTol, N,      &
                    B, HX, PX, X, ERR_RES, id)      
      !if (id==rootID) print *, ' QMR Iter:', OSB_QMR
  end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMR_OSB( N, B, X, ERR_RES)
      integer, intent(IN)     :: N
      double precision, intent(in)  :: B(N) 
      double precision, intent(out) :: X(N) 
      double precision, intent(out) :: ERR_RES

      integer :: qmr_mpi

      QMR_OSB = QMR_MPI(sQMRConvType, sQMR%mMax, sQMR%mTol, N,      &
                    B, EHX, EPX, X, ERR_RES, id)

  end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMR_OSB0( N, B, X, ERR_RES)
      integer, intent(IN)     :: N
      double precision, intent(in)  :: B(N) 
      double precision, intent(out) :: X(N) 
      double precision, intent(out) :: ERR_RES

      integer :: qmr_mpi

      QMR_OSB0 = QMR_MPI(sQMRConvType, sQMR%mMax, sQMR%mTol, N,      &
                    B, EHX, EPX0, X, ERR_RES, id)

  end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMR_OSBD1( N, B, X, ERR_RES)
      integer, intent(IN)     :: N
      double precision, intent(in)  :: B(N) 
      double precision, intent(out) :: X(N) 
      double precision, intent(out) :: ERR_RES

      integer :: qmr_MPI

      QMR_OSBD1 = QMR_MPI(sQMRConvType, sQMR%mMax, sQMR%mTol, N,   &
                     B, EHX, EPXD1, X, ERR_RES, id)

  end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMR_OSBD2( N, B, X, ERR_RES)
      integer, intent(IN)     :: N
      double precision, intent(in)  :: B(N) 
      double precision, intent(out) :: X(N) 
      double precision, intent(out) :: ERR_RES

      integer :: qmr_mpi

      QMR_OSBD2 = QMR_mpi(sQMRCONVTYPE, sQMR%mMAX, sQMR%mTol, N,    &
                     B, EHX, EPXD2, X, ERR_RES, id)

  end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    The following functions don't allocate memory for OSBW     c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMR_OSBW( N, B,  X, ERES)  
      integer, intent(IN)     :: N
      double precision, intent(in) :: B(N) 
      double precision, intent(out):: X(N)
      double precision, intent(out):: ERES
      
      integer :: qmr_mpi
  
      QMR_OSBW = QMR_MPI(sQMRCONVTYPE,sQMR%mMAX,sQMR%mTOL,N,    &
               B, EHX, EPXW, X, ERES, id)

  end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c               Perform QMR using the first order of Hij              c
!c All the first order of Hij within the energy window are calculated  c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function QMR_OSBW_NOMEM( N, B,  X, ERES)  
      integer, intent(IN)     :: N
      double precision, intent(in)  :: B(N)
      double precision, intent(out) :: X(N)
      double precision, intent(out) :: ERES
      
      integer :: qmr_mpi
  
      if ( initMOSBW()) then
          QMR_OSBW_NOMEM = QMR_MPI(sQMRCONVTYPE,sQMR%mMAX,sQMR%mTOL, &
                         N, B, EHX, EPXW, X, ERES, id)
      else
          QMR_OSBW_NOMEM = 0
      end if
 
      call finalMOSBW()

  end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
