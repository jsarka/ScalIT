integer function PIST1_CX_MPI(MYID, ROOTID, E0, ETOL, N, X,      &
               M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  integer, intent(IN)  ::  MYID, ROOTID
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, M, M_MAX
  double complex, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv    
  double complex, intent(out) ::EIG(M)
  double precision, intent(out) :: RES               

  integer :: startcnt, stepcnt, stepeig, ntype
  integer :: PIST_CONV_CX_MPI

  startCnt = 2*M;  stepcnt  = M
  stepEIG  = 1;    ntype    = 0

  PIST1_CX_MPI = PIST_CONV_CX_MPI(MYID,ROOTID,E0,ETOL,ntype,N,X,        &
          startcnt,stepcnt,stepeig, M,M_MAX,H0X,LinSolv,EIG,RES)
  
end function 

!**********************************************************************************

integer function PIST2_CX_MPI(MYID, ROOTID, E0, ETOL, N, X,   &
                 M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  integer, intent(IN)  ::  MYID, ROOTID
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, M, M_MAX
  double complex, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv    
  double complex, intent(out) ::EIG(M)
  double precision, intent(out) :: RES               

  integer :: startcnt, stepcnt, stepeig, ntype
  integer :: PIST_CONV_CX_MPI

  startCnt = 2*M;  stepcnt  = M
  stepEIG  = M;     ntype    = 0

  PIST2_CX_MPI=PIST_CONV_CX_MPI(MYID,ROOTID,E0,ETOL,ntype,N,X,   &
                 startcnt, stepcnt, stepeig, M,M_MAX,H0X, LinSolv,EIG,RES)
  
end function 
!**********************************************************************************


!**********************************************************************************
integer function PIST51_CX_MPI( MYID, ROOTID,E0, ETOL, N, X,   &
                 M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  integer, intent(IN)  ::  MYID, ROOTID
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, M, M_MAX
  double complex, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv
  double complex, intent(out) ::EIG(M)
  double precision, intent(out) :: RES               

  integer :: startcnt, stepcnt, stepeig, ntype
  integer :: PIST_CONVERG_CX_MPI

  startCnt = 2*M;  stepcnt  = M
  stepEIG  = 1;    ntype    = 0

  PIST51_CX_MPI = PIST_CONVERG_CX_MPI( MYID, ROOTID,E0,ETOL,ntype, N,X, &
          startcnt,stepcnt, stepeig, M,M_MAX,H0X, LinSolv,EIG,RES)
  
end function 

!**********************************************************************************

!**********************************************************************************
integer function PIST52_CX_MPI( MYID, ROOTID,E0, ETOL, N, X,    &
                M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  integer, intent(IN)  ::  MYID, ROOTID
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, M, M_MAX
  double complex, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv
  double complex,  intent(out) ::EIG(M)  
  double precision, intent(out) :: RES               

  integer :: startcnt, stepcnt, stepeig, ntype
  integer :: PIST_CONVERG_CX_MPI

  startCnt = 2*M;  stepcnt  = M
  stepEIG  = M;    ntype    = 0

  PIST52_CX_MPI=PIST_CONVERG_CX_MPI(MYID,ROOTID,E0,ETOL,ntype,N,X,    &
         startcnt, stepcnt, stepeig, M,M_MAX,H0X, LinSolv, EIG,RES)
  
end function 
!**********************************************************************************
