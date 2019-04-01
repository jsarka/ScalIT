integer function PIST1_SX(E0, ETOL, N, X, M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL   
  integer, intent(IN)  :: N, M, M_MAX
  double complex, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv
  double complex, intent(out) ::EIG(M)  
  double precision, intent(out) :: RES               

  integer :: startcnt, stepcnt, stepeig, ntype
  integer :: PIST_CONV_SX

  startCnt = 2*M;  stepcnt  = M
  stepEIG  = 1;    ntype    = 0

  PIST1_SX   = PIST_CONV_SX(E0,ETOL,ntype, N,X,startcnt,      &
                stepcnt, stepeig, M,M_MAX,H0X, LinSolv, EIG,RES)
  
end function 

!*****************************************************************************

integer function PIST2_SX(E0, ETOL, N, X, M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, M, M_MAX
  double complex, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv    
  double complex, intent(out) ::EIG(M)  
  double precision, intent(out) :: RES               

  integer :: startcnt, stepcnt, stepeig, ntype
  integer :: PIST_CONV_SX

  startCnt = 2*M;  stepcnt  = M
  stepEIG  = M;    ntype    = 0

  PIST2_SX   = PIST_CONV_SX(E0,ETOL,ntype, N,X,startcnt,        &
                stepcnt, stepeig, M,M_MAX,H0X, LinSolv,EIG,RES)
  
end function 
!*****************************************************************************


!*****************************************************************************
integer function PIST51_SX(E0, ETOL, N, X, M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, M, M_MAX
  double complex, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv
  double complex, intent(out) ::EIG(M)  
  double precision, intent(out) :: RES               

  integer :: startcnt, stepcnt, stepeig, ntype
  integer :: PIST_CONVERG_SX

  startCnt = 2*M;  stepcnt  = M
  stepEIG  = 1;    ntype    = 0

  PIST51_SX   = PIST_CONVERG_SX(E0,ETOL,ntype, N,X,startcnt,        &
                stepcnt, stepeig, M,M_MAX,H0X, LinSolv, EIG,RES)
  
end function 

!****************************************************************************

!*****************************************************************************
integer function PIST52_SX(E0, ETOL, N, X, M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL 
  integer, intent(IN)  :: N, M, M_MAX
  double complex, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv    
  double complex, intent(out) ::EIG(M)
  double precision, intent(out) :: RES               

  integer :: startcnt, stepcnt, stepeig, ntype
  integer :: PIST_CONVERG_SX

  startCnt = 2*M;  stepcnt  = M
  stepEIG  = M;    ntype    = 0

  PIST52_SX   = PIST_CONVERG_SX(E0,ETOL,ntype, N,X,startcnt,        &
                stepcnt, stepeig, M,M_MAX,H0X, LinSolv,EIG,RES)
  
end function 
!*****************************************************************************
