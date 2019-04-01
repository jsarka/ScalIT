
integer function PIST1_MPI(MYID, ROOTID, E0, ETOL, N, X,   &
               M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  integer, intent(IN)  ::  MYID, ROOTID
  double precision, intent(in) :: E0, ETOL   
  integer, intent(IN)  :: N, M, M_MAX
  double precision, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv    
  double precision, intent(out) ::EIG(M)  
  double precision, intent(out) :: RES               

  integer :: startcnt, stepcnt, stepeig, ntype
  integer :: PIST_CONV_MPI

  startCnt = 2*M;  stepcnt  = M
  stepEIG  = 1;    ntype    = 0

  PIST1_MPI=PIST_CONV_MPI(MYID, ROOTID, E0, ETOL, ntype, N, X,          &
                 startcnt,stepcnt,stepeig, M,M_MAX,H0X,LinSolv,EIG,RES)
  
end function 

!**********************************************************************************

integer function PIST2_MPI( MYID, ROOTID, E0, ETOL, N, X,   &
                 M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  integer, intent(IN)  ::  MYID, ROOTID
  double precision, intent(in) :: E0, ETOL 
  integer, intent(IN)  :: N, M, M_MAX
  double precision, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv
  double precision, intent(out) ::EIG(M)  
  double precision, intent(out) :: RES               

  integer :: startcnt, stepcnt, stepeig, ntype
  integer :: PIST_CONV_MPI

  startCnt = 2*M;  stepcnt  = M
  stepEIG  = M;    ntype    = 0

  PIST2_MPI = PIST_CONV_MPI(MYID, ROOTID,E0,ETOL,ntype, N,X, &
                 startcnt, stepcnt, stepeig, M,M_MAX,H0X, LinSolv,EIG,RES)
  
end function 
!**********************************************************************************


!**********************************************************************************
integer function PIST51_MPI(MYID, ROOTID,E0, ETOL, N, X,  &
                 M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  integer, intent(IN)  ::  MYID, ROOTID
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, M, M_MAX
  double precision, intent(in) :: X (N)
  external             :: H0X
  integer, external    :: LinSolv    
  double precision, intent(out) ::EIG(M)
  double precision, intent(out) :: RES               

  integer :: startcnt, stepcnt, stepeig, ntype
  integer :: PIST_CONVERG_MPI

  startCnt = 2*M;   stepcnt  = M
  stepEIG  = 1;     ntype    = 0

  PIST51_MPI = PIST_CONVERG_MPI( MYID, ROOTID,E0,ETOL,ntype, N,X,    &
          startcnt,stepcnt, stepeig, M,M_MAX,H0X,LinSolv,EIG,RES)
  
end function 

!**********************************************************************************

!**********************************************************************************
integer function PIST52_MPI(MYID, ROOTID,E0, ETOL, N, X,   &
                M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  integer, intent(IN)  :: MYID, ROOTID
  double precision, intent(in) :: E0, ETOL 
  integer, intent(IN)  :: N, M, M_MAX
  double precision, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv
  double precision, intent(out) ::EIG(M)  
  double precision, intent(out) :: RES               

  integer :: startcnt, stepcnt, stepeig, ntype
  integer :: PIST_CONVERG_MPI

  startCnt = 2*M;  stepcnt  = M
  stepEIG  = M;    ntype    = 0

  PIST52_MPI = PIST_CONVERG_MPI(MYID, ROOTID,E0,ETOL,ntype, N,X, &
               startcnt, stepcnt, stepeig, M,M_MAX,H0X, LinSolv,EIG,RES)
  
end function 
!**********************************************************************************
