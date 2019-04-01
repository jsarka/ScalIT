!*****************************************************************
!*     COMPLEX VERSION OF PIST IN SEQUENTIAL ENVIRONMENT         *
!*****************************************************************

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A fast version but difficult to implement                  c
!c         It will reuse the previous generated VJ, and Hij               c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function PISTF_CONV_CX(E0, ETOL, NTYPE, N, X, startCNT, stepCnt,  &
                  stepeig, M, M_MAX, H0XCX, LinSolvCX, EIG, RES, Filename)
  implicit none
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, NTYPE
  double complex, intent(in) :: X(N) 
  integer, intent(IN)  :: STARTCNT, STEPCNT, STEPEIG        
  integer, intent(IN)  :: M, M_MAX 
  external             :: H0XCX
  integer, external    :: LinSolvCX     
  character(LEN=*)     :: filename

  double complex, intent(out) ::EIG(M)   
  double precision, intent(out) :: RES              

  integer, parameter  :: nDiffTYPE    = 1  ! Real part
  integer, parameter  :: nConvType    = 6  ! Real^2+Img^2

!cccccccccccccccccccccccccccccccc
  double complex :: VJ(N, M_MAX),HMat(M_Max,M_Max),Hmat_old(M_MAX,M_MAX) 
  double complex :: tmpEig(M_MAX),oldEig(M)  
  
  integer   :: cnt_high, cnt_low  
  integer   :: I, J, NUM

!cccccccccccccccccccccccccccccccc  
  double precision  :: getDiffMax_CX, getSepMax_cx
  integer           :: PISTH0_CX

!ccccccccccccc   Parameters to call ZGEEV subroutine cccccccccccccc
  integer, parameter  :: LDVL=1,  LDVR=1, SCALE = 3
  double complex      :: VL(LDVL), VR(LDVR), WORK(SCALE*M_Max), RWORK(2*N)
  integer             ::  LWORK 
 
  LWORK = SCALE*M_MAX  

!***************************************************** 

  if ((stepEig >= M_MAX) .or. (startCnt > M_MAX) .or. (M > M_MAX)) then
      PISTF_CONV_CX = 0;      return
  end if

  cnt_low = 1;   cnt_high  = max(startCnt, M, stepEig+1) 

  do  while ( cnt_high <=M_MAX )

     NUM=PISTH0_CX(N,X,M_MAX,CNT_LOW,CNT_HIGH,H0XCX,LinSolvCX,VJ,HMAT_OLD)

     if ( NUM <= 0) then
         PISTF_CONV_CX = -CNT_HIGH;    return
     end if  

     num = cnt_high-stepEig
     HMAT(1:num,1:num) = HMAT_OLD(1:num,1:num)  
     call ZGEEV('N', 'N', NUM, HMAT, M_MAX, tmpEig, VL, LDVL, VR, LDVR,    &
                 WORK, LWORK, RWORK, PISTF_Conv_CX)
     if (PISTF_CONV_CX /= 0) then       ! Lapack error   
         PISTF_CONV_CX = -NUM;        return                
     end if
     call getWindow_CX(E0, nDiffType, NUM, tmpEig, M, oldEig) 

     NUM = cnt_high
     HMAT(1:NUM,1:NUM) = HMAT_old(1:NUM,1:NUM)
     call ZGEEV('N', 'N', NUM, HMAT, M_MAX, tmpEIG,VL, LDVL, VR, LDVR,       &
                WORK, LWORK, RWORK, PISTF_Conv_CX)
     if (PISTF_CONV_CX /= 0) then       ! Lapack error   
         PISTF_CONV_CX = -NUM;      return                
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
          PISTF_CONV_CX = cnt_high;      exit
     end if    

     PISTF_CONV_CX = - cnt_high 

     cnt_low  = cnt_high + 1
     if (cnt_low >= M_MAX) then
           pistf_conv_CX = - M_MAX;        exit
      end if

      cnt_high = cnt_high + stepcnt
      if (cnt_high > M_MAX) cnt_high = M_MAX
  end do

  if (PISTF_CONV_CX>0) then  ! store data
     open(99, FILE=filename, status='Replace', FORM='UNFORMATTED')
     write(99) PISTF_CONV_CX, N
     write(99) tmpEig(1:PISTF_CONV_CX),HMat_old(1:PISTF_CONV_CX,1:PISTF_CONV_CX)
     write(99) VJ(1:N, 1:PISTF_CONV_CX)
     close(99)
  end if

  call Reorder_CX('A', nDiffType,  M, eig)
  
end function PISTF_CONV_CX
!**************************************************************************

