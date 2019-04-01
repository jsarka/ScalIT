integer function PIST1_DX(E0, ETOL, N, X, M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL   
  integer, intent(IN)  :: N, M, M_MAX
  double precision, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv
  double complex, intent(out) ::EIG(M)  
  double precision, intent(out) :: RES               

  integer :: startcnt, stepcnt, stepeig, ntype
  integer :: PIST_CONV_DX

  startCnt = 2*M;  stepcnt  = M
  stepEIG  = 1;    ntype    = 0

  PIST1_DX   = PIST_CONV_DX(E0,ETOL,ntype, N,X,startcnt,      &
                stepcnt, stepeig, M,M_MAX,H0X, LinSolv, EIG,RES)
  
end function 

!*****************************************************************************

integer function PIST2_DX(E0, ETOL, N, X, M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, M, M_MAX
  double precision, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv    
  double complex, intent(out) ::EIG(M)  
  double precision, intent(out) :: RES               

  integer :: startcnt, stepcnt, stepeig, ntype
  integer :: PIST_CONV_DX

  startCnt = 2*M;  stepcnt  = M
  stepEIG  = M;    ntype    = 0

  PIST2_DX   = PIST_CONV_DX(E0,ETOL,ntype, N,X,startcnt,        &
                stepcnt, stepeig, M,M_MAX,H0X, LinSolv,EIG,RES)
  
end function 
!*****************************************************************************


!*****************************************************************************
integer function PIST51_DX(E0, ETOL, N, X, M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, M, M_MAX
  double precision, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv
  double complex, intent(out) ::EIG(M)  
  double precision, intent(out) :: RES               

  integer :: startcnt, stepcnt, stepeig, ntype
  integer :: PIST_CONVERG_DX

  startCnt = 2*M;  stepcnt  = M
  stepEIG  = 1;    ntype    = 0

  PIST51_DX   = PIST_CONVERG_DX(E0,ETOL,ntype, N,X,startcnt,        &
                stepcnt, stepeig, M,M_MAX,H0X, LinSolv, EIG,RES)
  
end function 

!****************************************************************************

!*****************************************************************************
integer function PIST52_DX(E0, ETOL, N, X, M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL 
  integer, intent(IN)  :: N, M, M_MAX
  double precision, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv    
  double complex, intent(out) ::EIG(M)
  double precision, intent(out) :: RES               

  integer :: startcnt, stepcnt, stepeig, ntype
  integer :: PIST_CONVERG_DX

  startCnt = 2*M;  stepcnt  = M
  stepEIG  = M;    ntype    = 0

  PIST52_DX   = PIST_CONVERG_DX(E0,ETOL,ntype, N,X,startcnt,        &
                stepcnt, stepeig, M,M_MAX,H0X, LinSolv,EIG,RES)
  
end function 
!*****************************************************************************
