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
!c   double(N)  X  : intial vector                                   c
!c   integer M     : number of eigen values                          c
!c   integer M_MAX : max. number of Lanczos iterations               c
!c   external HX   : function for H*X(N, v, w) , w=H*v               c
!c                                                                   c
!c Output parameters:                                                c
!c   double(M) eig : final eigen values                              c
!c   double    res : max. error for the convergence testing          c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function LANCONV1(E0, ETOL, N, X, M, M_MAX, HX, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, M, M_MAX               
  double precision, intent(in) :: X(N)
  external             :: HX                       
  double precision, intent(out) ::EIG(M)   
  double precision, intent(out) :: RES                

  integer :: startCnt, stepCnt, stepEig, nType
  integer :: LAN_CONV

  startCnt = max(int(DSQRT(1.0D0*M_MAX)), M)
  stepCnt  = startCnt;   stepEig  = 1;   nType    = 0

  LANCONV1   = LAN_CONV(E0, ETOL, nType, N, X, startCnt, stepCnt, stepEig,   &
                        M, M_MAX, HX, EIG, RES)
  
end function LANCONV1

!*****************************************************************************
integer function LANCONV2(E0, ETOL, N, X, M, M_MAX, HX, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, M, M_MAX               
  double precision, intent(in) :: X(N)
  external             :: HX                       
  double precision, intent(out) ::EIG(M)   
  double precision, intent(out) :: RES                

  integer :: startCnt, stepCnt, stepEig, nType
  integer :: LAN_CONV

  startCnt = max(int(DSQRT(1.0D0*M_MAX)), M)
  stepCnt  = startCnt;    stepEig  = startCnt/2
  nType    = 0

  LANCONV2   = LAN_CONV(E0, ETOL, nType, N, X, startCnt, stepCnt, stepEig,    &
                        M, M_MAX, HX, EIG, RES)
  
end function LANCONV2
!******************************************************************************

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                   c
!c   Using the default start and step sizes for the slow Lanczos     c
!c            It needs less memory to store the vector               c
!c                                                                   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function LANCONV51(E0, ETOL, N, X, M, M_MAX,HX, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, M, M_MAX               
  double precision, intent(in) :: X(N)
  external             :: HX                       
  double precision, intent(out) ::EIG(M)
  double precision, intent(out) :: RES                

  integer :: startCnt, stepCnt, stepEig, nType
  integer :: LAN_CONVERG

  startCnt = max(int(DSQRT(1.0D0*M_MAX)), M)
  stepCnt  = startCnt;    stepEig  = 1;    nType    = 0

  LANCONV51 = LAN_CONVERG(E0, ETOL, nType, N, X, startCnt, stepCnt, stepEig, &
                           M, M_MAX, HX, EIG, RES)
  
end function LANCONV51

!******************************************************************************
integer function LANCONV52(E0, ETOL, N, X, M, M_MAX, HX, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, M, M_MAX               
  double precision, intent(in) :: X(N)
  external             :: HX                       
  double precision, intent(out) ::EIG(M)   
  double precision, intent(out) :: RES                

  integer :: startCnt, stepCnt, stepEig, nType
  integer :: LAN_CONVERG

  startCnt = max(int(DSQRT(1.0D0*M_MAX)), M)
  stepCnt  = startCnt;   stepEig  = startCnt/2
  nType    = 0

  LANCONV52 = LAN_CONVERG(E0, ETOL, nTYPE, N, X, startCnt, stepCnt, stepEig, &
                          M, M_MAX, HX, EIG, RES)
  
end function LANCONV52
!******************************************************************************
