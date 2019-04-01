!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                   c
!c   Using the default start and step sizes for the fast Lanczos     c
!c            It needs more memory to store the vector               c
!c                                                                   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                   c
!c Input Parameters:                                                 c
!c   double  E0    : central energy                                  c
!c   double  Etol  : energy tolerence for convergence                c
!c   integer N     : dimension of initial vector                     c
!c   complex(N)  x : intial vector                                   c
!c   integer M     : number of eigen values                          c
!c   integer M_MAX : max. number of Lanczos iterations               c
!c   external HXCX : function for H*X(N, v, w) , w=H*v               c
!c                                                                   c
!c Output parameters:                                                c
!c   double(M) eig : final eigen values                              c
!c   double    res : max. error for the convergence testing          c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function LANCONV1_DX(E0, ETOL, N, X, M, M_MAX, HXCX, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL              
  integer, intent(IN)  :: N , M, M_MAX                     
  double complex,   intent(in) :: X(N)
  external             :: HXCX         
  double precision, intent(out) ::EIG(M)
  double precision, intent(out) :: RES               

  integer :: startCnt, stepCnt, stepEig, nType
  integer :: LAN_CONV_DX

  startCnt = max(int(DSQRT(1.0D0*M_MAX)), M)
  stepCnt  = startCnt;  stepEig  = 1;  nType    = 0

  LANCONV1_DX = LAN_CONV_DX( E0, ETOL, nType, N, X, startCnt, stepCnt,       &
                      stepEig, M, M_MAX, HXCX, EIG, RES)
  
end 

!****************************************************************************
integer function LANCONV2_DX(E0, ETOL, N, X, M, M_MAX, HXCX, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL              
  integer, intent(IN)  :: N , M, M_MAX                     
  double complex,  intent(in) :: X(N)
  external             :: HXCX         
  double precision, intent(out) ::EIG(M)
  double precision, intent(out) :: RES    

  integer :: startCnt, stepCnt, stepEig, nType
  integer :: LAN_CONV_DX

  startCnt = max(int(DSQRT(1.0D0*M_MAX)), M)
  stepCnt  = startCnt;   stepEig = startCnt/2;   nType = 0

  LANCONV2_DX = LAN_CONV_DX(E0, ETOL, nType, N, X, startCnt, stepCnt,     &
                            stepEig, M, M_MAX, HXCX, EIG, RES)
  
end
!**************************************************************************

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                   c
!c   Using the default start and step sizes for the slow Lanczos     c
!c            It needs less memory to store the vector               c
!c                                                                   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function LANCONV51_DX(E0, ETOL, N, X, M, M_MAX,HXcx, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL              
  integer, intent(IN)  :: N , M, M_MAX                     
  double complex,   intent(in) :: X(N)
  external             :: HXCX         
  double precision, intent(out) ::EIG(M)
  double precision, intent(out) :: RES    

  integer :: startCnt, stepCnt, stepEig, nType
  integer :: LAN_CONVERG_DX

  startCnt = max(int(DSQRT(1.0D0*M_MAX)), M)
  stepCnt  = startCnt;   stepEig = 1;    nType = 0

  LANCONV51_DX = LAN_CONVERG_DX(E0, ETOL, nType, N, X, startCnt, stepCnt,   &
                                stepEig, M, M_MAX, HXCX, EIG, RES)
  
end 

!****************************************************************************
integer function LANCONV52_DX(E0, ETOL, N, X, M, M_MAX, HXcx, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL              
  integer, intent(IN)  :: N , M, M_MAX                     
  double complex,  intent(in) :: X(N)  
  external             :: HXCX         
  double precision, intent(out) ::EIG(M)
  double precision, intent(out) :: RES    

  integer :: startCnt, stepCnt, stepEig, nType
  integer :: LAN_CONVERG_DX

  startCnt = max(int(DSQRT(1.0D0*M_MAX)), M)
  stepCnt  = startCnt;   stepEig = startCnt/2;   nType = 0

  LANCONV52_DX = LAN_CONVERG_DX(E0, ETOL, nTYPE, N, X, startCnt, stepCnt,   &
                                stepEig, M, M_MAX, HXCX, EIG, RES)
  
end 
!****************************************************************************
