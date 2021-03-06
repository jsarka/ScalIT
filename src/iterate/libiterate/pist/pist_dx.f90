!*****************************************************************
!*     COMPLEX VERSION OF PIST IN SEQUENTIAL ENVIRONMENT         *
!*****************************************************************
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A fast version but difficult to implement                  c
!c         It will reuse the previous generated VJ, and Hij               c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function PIST_CONV_DX(E0, ETOL, NTYPE, N, X, startCNT,   &
                  STEPcnt, stepeig, M, M_MAX, H0XDX, LinSolv, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, NTYPE
  double precision, intent(in) :: X(N) 
  integer, intent(IN)  :: STARTCNT, STEPCNT, STEPEIG        
  integer, intent(IN)  :: M, M_MAX 
  external             :: H0XDX
  integer, external    :: LinSolv
  double complex,  intent(out) ::EIG(M)   
  double precision,intent(out) :: RES              

  integer, parameter   :: nDiffTYPE    = 1  ! Real part
  integer, parameter   :: nConvType    = 6  ! Real^2+Img^2

!cccccccccccccccccccccccccccccccc
  double precision :: VJ(N, M_max)
  double complex :: Hmat(M_max,M_max),HMAT_OLD(M_max,M_max)
  double complex :: tmpEig(M_MAX), oldEig(M)  
  
  integer   :: cnt_high, cnt_low  
  integer   :: I, J, NUM

!cccccccccccccccccccccccccccccccc  
  double precision  :: getDiffMax_CX, getSepMax_cx
  integer           :: PISTH0_DX

!ccccccccccccc   Parameters to call ZGEEV subroutine cccccccccccccc
  integer, parameter   ::  LDVL=1, LDVR=1,SCALE=3
  double complex  :: VL(LDVL), VR(LDVR),WORK(SCALE*M_max),RWork(2*N) 
  integer         :: LWORK 
  
  LWORK = SCALE*M_MAX  

!************************************ 
  if ((stepEig >= M_MAX) .or. (startCnt > M_MAX) .or. (M > M_MAX)) then
      PIST_CONV_DX = 0;     return
  end if

  cnt_low = 1;   cnt_high = max(startCnt, M, stepEig+1) 

  do  while ( cnt_high <=M_MAX )

     NUM=PISTH0_DX(N,X,M_MAX,CNT_LOW,CNT_HIGH,H0XDX,LinSolv,VJ,HMAT_OLD)

     if ( NUM <= 0) then
         PIST_CONV_DX = -CNT_HIGH;   return
     end if  

     num = cnt_high-stepEig
     HMAT(1:num,1:num) = HMAT_OLD(1:num,1:num)  
     call ZGEEV('N', 'N', NUM, HMAT, M_MAX, tmpEig, VL, LDVL, VR, LDVR,    &
                 WORK, LWORK, RWORK, PIST_Conv_DX)
     if (PIST_CONV_DX /= 0) then       ! Lapack error   
         PIST_CONV_DX = -NUM;    return                
     end if
     call getWindow_CX(E0, nDiffType, NUM, tmpEig, M, oldEig) 

     num = cnt_high
     HMAT(1:NUM,1:NUM) = HMAT_old(1:NUM,1:NUM)
     call ZGEEV('N', 'N', NUM, HMAT, M_MAX, tmpEIG,VL, LDVL, VR, LDVR,       &
                WORK, LWORK, RWORK, PIST_Conv_DX)
     if (PIST_CONV_DX /= 0) then       ! Lapack error   
         PIST_CONV_DX = -NUM;    return                
     end if

     call getWindow_CX(E0, nDiffType, NUM, tmpEig, M, Eig) 

     select case (nType)
     case (:0)
          RES = getDiffMax_CX(nConvTYPE, M,eig, oldEig)        

     case (1:)
          RES = getSepMax_CX(nConvTYPE, M, oldEig, cnt_high, tmpEig)
     end select

     print *, 'Cnt1=',CNT_HIGH, ' Cnt2=',(cnt_High-stepEig),' RES=', RES

     if ( RES <= ETOL ) then
          PIST_CONV_DX = cnt_high;        exit
     end if    

     PIST_CONV_DX = - cnt_high 

     cnt_low  = cnt_high + 1
     if (cnt_low >= M_MAX) then
           pist_conv_DX = - M_MAX;      exit
      end if

      cnt_high = cnt_high + stepcnt
      if (cnt_high > M_MAX)  cnt_high = M_MAX

  end do

  call Reorder_CX('A', nDiffType,  M, eig)
  
end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A slow version but easy to implement                       c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function PIST_CONVERG_DX(E0, ETOL, NTYPE, N, X, startCNT,   &
                  STEPcnt, stepeig, M, M_MAX, H0XDX, LinSolv, EIG, RES)
  implicit none
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, NTYPE
  double precision, intent(in) :: X(N) 
  integer, intent(IN)  :: STARTCNT, STEPCNT, STEPEIG        
  integer, intent(IN)  :: M, M_MAX 
  external             :: H0XDX
  integer, external    :: LinSolv     

  double complex,  intent(out) ::EIG(M)   
  double precision, intent(out) :: RES              

  integer, parameter   :: nDiffTYPE    = 1  ! Real part
  integer, parameter   :: nConvType    = 6  ! Real^2+Img^2

!cccccccccccccccccccccccccccccccc
  double complex :: tmp1Eig(M_MAX), tmp2Eig(M_MAX), oldEig(M)
  integer        ::  cnt1, cnt2

!cccccccccccccccccccccccccccccccc
  integer           :: PIST_EIG_DX
  double precision  :: getDiffMax_CX, getSepMax_CX
   
!************************************  
  if ((stepEig >= M_MAX) .or. (startCnt > M_MAX) .or. (M > M_MAX))  then
      PIST_CONVERG_DX = 0;      return
  end if

  Cnt1  = max(startCnt, M, stepEig+1) 

  PIST_CONVERG_DX = - M_MAX

  do  while ( CNT1 <=M_MAX )
      CNT2    = CNT1-STEPEIG

      PIST_CONVERG_DX=PIST_EIG_DX(N,X,CNT1,CNT2,H0XDX,LinSolv,tmp1EIG,tmp2Eig)

      if (PIST_CONVERG_DX <= 0) then
          PIST_CONVERG_DX = -CNT1;        return
      end if

      call getWindow_CX(E0, nDiffType, CNT1, tmp1Eig, M, Eig) 
      call getWindow_CX(E0, nDiffType, CNT2, tmp2Eig, M, oldEig)     
      
      select case (nType)
      case (:0)
          RES = getDiffMax_CX(nConvType,M, eig, oldEig)        

      case (1:)
           RES = getSepMax_CX(nConvTYPE, M, oldEig, cnt1, tmp1Eig)
      end select  

      print *, 'Cnt1=', CNT1,' Cnt2=', CNT2, ' RES=', RES
  
      if ( RES < ETOL ) then
          PIST_CONVERG_DX = CNT1;     exit
      end if    

      CNT1 = CNT1 + stepCnt

      if (cnt1>M_MAX) cnt1=M_MAX
  end do

  call Reorder_CX('A', nDIFFTYPE, M, eig)
  
end function PIST_CONVERG_DX
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    PIST algorithm  to get the eigen values near energy E        c
!c       Apply for Symmetric matrix H                              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function PIST_DX(N, X, M, H0XDX, LinSolv, EIG)
  implicit none
  integer, intent(IN)  :: N,M
  double precision, intent(in) :: X(N)  
  external             :: H0XDX
  integer, external    :: LinSolv
  double complex, intent(out) :: EIG(M)

!ccccccccccccccc  Local Variables cccccccccccccccccccccccccc  
  double complex :: HMAT(M,M)   ! HMAT = <U|H|U>
  integer :: pistHij_DX

  !ccccccccccccc   Parameters to call ZGEEV subroutine cccccccccccccc
  integer, parameter   :: LDVL=1, LDVR=1, SCALE=3
  double complex       :: VL(LDVL), VR(LDVR), Work(SCALE*M), RWork(2*N)
  integer              :: LWORK 

  LWORK = SCALE * M  

  PIST_DX = PISTHIJ_DX(N, X, M, H0XDX, LinSolv, HMAT)

  if (PIST_DX > 0) then
      call ZGEEV('N', 'N', M, HMAT, M, EIG, VL, LDVL, VR, LDVR, WORK,   &
                   LWORK, RWORK, PIST_DX)
      if (PIST_DX == 0) then      
          PIST_DX = M
      else
          PIST_DX = -M    
      end if
  end if

end 
!************************************************************************


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    PIST algorithm  to get the eigen values near energy E        c
!c      For two sets of,  Apply for Symmetric matrix H             c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function PIST_EIG_DX(N, X, M1, M2, H0XDX, LinSolv, EIG1, EIG2)
  implicit none
  integer, intent(IN)  :: N, M1, M2
  double precision, intent(in) :: X(N)  
  external             :: H0XDX
  integer, external    :: LinSolv
  double complex, intent(out) :: EIG1(M1), EIG2(M2)

!ccccccccccccccc  Local Variables cccccccccccccccccccccccccc  
  integer  :: M_MAX, INFO, LWORK
  double complex, allocatable :: HMAT(:,:), HOLD(:,:)   ! HMAT = <U|H|U>
  double complex, allocatable :: WORK(:), RWORK(:)

  integer, parameter ::  LDVL = 1, LDVR=1, SCALE=3
  double complex  :: VL(LDVL), VR(LDVR) 
  integer      ::  PISTHIJ_DX
  
  M_MAX = max(M1, M2);   LWORK = SCALE*M_MAX  
  PIST_EIG_DX = -1

  allocate( HMAT(M_MAX, M_MAX),HOLD(M_MAX, M_MAX),        &
            WORK(LWORK),RWORK(2*M_MAX), STAT=INFO)

  if (INFO /= 0)   return 
 
  PIST_EIG_DX = PISTHIJ_DX(N, X, M_MAX, H0XDX, LinSolv, HMAT)

  if ( PIST_EIG_DX > 0 ) then
      HOLD(1:M_MAX, 1:M_MAX) = HMAT(1:M_MAX, 1:M_MAX)
      call ZGEEV('N', 'N', M1, HMAT, M_MAX, EIG1, VL, LDVL, VR, LDVR, WORK, &
                 LWORK, RWORK, INFO)

      if (INFO /= 0) then
           PIST_EIG_DX = - M_MAX
      else
           call ZGEEV('N', 'N', M2, HOLD, M_MAX, EIG2, VL, LDVL, VR, LDVR,  &
                 WORK, LWORK, RWORK, INFO)

           if (INFO /= 0)   PIST_EIG_DX = - M_MAX
      end if
  end if

  deallocate(HOLD, HMAT, WORK, RWORK,STAT=INFO)

end 
!**************************************************************************

