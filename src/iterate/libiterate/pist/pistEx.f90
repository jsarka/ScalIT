
integer function PIST1(E0, ETOL, N, X, M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, M, M_MAX
  double precision, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv    
  double precision, intent(out) ::EIG(M)  
  double precision, intent(out) :: RES               

  integer :: startcnt, stepcnt, stepeig, ntype
  integer :: PIST_CONV

  startCnt = 2*M;   stepcnt  = M
  stepEIG  = 1;     ntype    = 0

  PIST1   = PIST_CONV(E0,ETOL,ntype, N,X,startcnt, stepcnt, stepeig,    &
                      M,M_MAX,H0X, LinSolv, EIG,RES)
  
end function 

!*****************************************************************************

integer function PIST2(E0, ETOL, N, X, M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, M, M_MAX
  double precision, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv
  double precision, intent(out) ::EIG(M)  
  double precision, intent(out) :: RES               

  integer :: startcnt, stepcnt, stepeig, ntype
  integer :: PIST_CONV

  startCnt = 2*M;   stepcnt  = M
  stepEIG  = M;     ntype    = 0

  PIST2   = PIST_CONV(E0,ETOL,ntype,N,X,startcnt,stepcnt,stepeig,M,     &
                  M_MAX,H0X,LinSolv,EIG,RES)
  
end function 
!****************************************************************************


!*****************************************************************************
integer function PIST51(E0, ETOL, N, X, M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, M, M_MAX
  double precision, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv    
  double precision, intent(out) ::EIG(M)  
  double precision, intent(out) :: RES               

  integer :: startcnt, stepcnt, stepeig, ntype
  integer :: PIST_CONVERG

  startCnt = 2*M;   stepcnt  = M
  stepEIG  = 1;     ntype    = 0

  PIST51   = PIST_CONVERG(E0,ETOL,ntype, N,X,startcnt,stepcnt, stepeig,     &
                          M,M_MAX,H0X, LinSolv,EIG,RES)
  
end function 

!*****************************************************************************

!*****************************************************************************
integer function PIST52(E0, ETOL, N, X, M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, M, M_MAX
  double precision, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv    
  double precision, intent(out) ::EIG(M)
  double precision, intent(out) :: RES               

  integer :: startcnt, stepcnt, stepeig, ntype
  integer :: PIST_CONVERG

  startCnt = 2*M;  stepcnt  = M
  stepEIG  = M;    ntype    = 0

  PIST52   = PIST_CONVERG(E0,ETOL,ntype, N,X,startcnt,stepcnt, stepeig,  &
                      M,M_MAX,H0X, LinSolv,EIG,RES)
  
end function 
!****************************************************************************
