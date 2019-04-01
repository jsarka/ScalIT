!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Using the default start and step sizes for the fast Lanczos     c
!c            It needs more memory to store the vector               c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function LANCONV1_MPI(MYID, ROOTID, E0, ETOL, N, X, M,  &
                              M_MAX, HX, EIG, RES)
  implicit none
  integer, intent(IN)   ::  MYID, ROOTID
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, M, M_MAX               
  double precision, intent(in) :: X(N)
  external             :: HX                       
  double precision, intent(out) ::EIG(M)   
  double precision, intent(out) :: RES                

  integer :: startCnt, stepCnt, stepEig, nType
  integer :: LAN_CONV_MPI

  startCnt = max(int(DSQRT(1.0D0*M_MAX)), M)
  stepCnt  = startCnt;  stepEig  = 1;   nType    = 0

  LANCONV1_MPI = LAN_CONV_MPI( MYID, ROOTID, E0, ETOL, nType, N, X, &
                     startCnt, stepCnt, stepEig, M, M_MAX, HX, EIG, RES)
  
end 

!*****************************************************************************
integer function LANCONV2_MPI( MYID, ROOTID, E0, ETOL, N, X, M,  &
                             M_MAX, HX, EIG, RES)
  implicit none
  integer, intent(IN)   ::  MYID, ROOTID
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, M, M_MAX               
  double precision, intent(in) :: X(N)
  external             :: HX                       
  double precision, intent(out) ::EIG(M)
  double precision, intent(out) :: RES                

  integer :: startCnt, stepCnt, stepEig, nType
  integer :: LAN_CONV_MPI

  startCnt = max(int(DSQRT(1.0D0*M_MAX)), M)
  stepCnt  = startCnt;  stepEig = startCnt/2;   nType = 0

  LANCONV2_MPI = LAN_CONV_MPI( MYID, ROOTID, E0, ETOL, nType, N, X, &
                     startCnt, stepCnt, stepEig, M, M_MAX, HX, EIG, RES)
  
end 
!******************************************************************************

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Using the default start and step sizes for the slow Lanczos     c
!c            It needs less memory to store the vector               c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function LANCONV51_MPI( MYID, ROOTID,E0, ETOL, N, X, M, &
                               M_MAX,HX, EIG, RES)
  implicit none
  integer, intent(IN)   ::  MYID, ROOTID
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, M, M_MAX               
  double precision, intent(in) :: X(N)
  external             :: HX                       
  double precision, intent(out) ::EIG(M)   
  double precision, intent(out) :: RES                

  integer :: startCnt, stepCnt, stepEig, nType
  integer :: LAN_CONVERG_MPI

  startCnt = max(int(DSQRT(1.0D0*M_MAX)), M)
  stepCnt  = startCnt;  stepEig = 1;  nType = 0

  LANCONV51_MPI = LAN_CONVERG_MPI( MYID, ROOTID, E0, ETOL, nType, N, &
                        X, startCnt, stepCnt, stepEig, M, M_MAX, HX, EIG, RES)
  
end 

!******************************************************************************
integer function LANCONV52_MPI( MYID, ROOTID, E0, ETOL, N, X, M,    &
                               M_MAX, HX, EIG, RES)
  implicit none
  integer, intent(IN)   ::  MYID, ROOTID
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, M, M_MAX               
  double precision, intent(in) :: X(N)
  external             :: HX                       
  double precision, intent(out) ::EIG(M)   
  double precision, intent(out) :: RES                

  integer :: startCnt, stepCnt, stepEig, nType
  integer :: LAN_CONVERG_MPI

  startCnt = max(int(DSQRT(1.0D0*M_MAX)), M)
  stepCnt  = startCnt;  stepEig = startCnt/2;  nType = 0

  LANCONV52_MPI = LAN_CONVERG_MPI( MYID, ROOTID, E0, ETOL, nTYPE, N, X,  &
                     startCnt, stepCnt, stepEig, M, M_MAX, HX, EIG, RES)
  
end 
!******************************************************************************
