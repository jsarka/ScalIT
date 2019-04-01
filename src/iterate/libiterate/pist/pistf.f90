!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A fast version but difficult to implement                  c
!c         It will reuse the previous generated VJ, and Hij               c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function PISTF_CONV(E0, ETOL, NTYPE, N, X, STARTCNT,  &
          STEPCNT, STEPEIG, M, M_MAX, H0X, LinSolv, EIG, RES, FileName)
  implicit none
  double precision, intent(in)  :: E0, ETOL
  integer, intent(IN )          :: nType, N
  double precision, intent(in)  :: X(N)  
  integer, intent(IN)           :: startCNT, stepCNT, stepEIG
  integer, intent(IN)           :: M, M_MAX
  external                      :: H0X
  integer, external             :: LinSolv   
  double precision, intent(out) :: EIG(M)   
  double precision, intent(out) :: RES                
  character(LEN=*), intent(IN)  :: filename

!cccccccccccccccccccccccccccccccc
  double precision :: VJ(N, M_MAX) 
  double precision :: HMat(M_MAX,M_MAX), HMat_Old(M_MAX, M_MAX)
  double precision :: tmpEig(M_MAX), oldEig(M)  
  
  integer :: cnt_high, cnt_low, cnt_start
  integer :: I, J, NUM

!cccccccccccccccccccccccccccccccc  
  double precision  :: getDiffMax, getSepMax
  integer           :: PISTH0

!ccccccccccccc   Parameters to call DSYEV subroutine cccccccccccccc
  integer   :: LWORK  
  double precision :: Work(3*M_MAX) 
  LWORK = 3*M_MAX  

!************************************ 
  if ((stepEig >= M_MAX) .OR. (startCnt > M_MAX) .OR. (M > M_MAX)) then
      PISTF_CONV = 0;     return
  end if

  cnt_low   = 1;     cnt_high  = max(startCnt, M, stepEig+1)   

  do  while ( cnt_high <=M_MAX )

     PISTF_CONV = PISTH0(N, X,M_MAX,CNT_LOW,CNT_HIGH,H0X, LinSolv, VJ,HMAT_OLD)

     if (PISTF_CONV < 0) then
         PISTF_CONV = - CNT_HIGH;     return
     end if    

     num = cnt_high-stepEig
     HMAT(1:num,1:num) = HMAT_OLD(1:num,1:num)     
     call DSYEV('N', 'U', num, HMAT, M_MAX, tmpEIG, WORK, LWORK, PISTF_Conv)

     if (PISTF_CONV /= 0) then       ! Lapack error   
         PISTF_CONV = -num;         return                
     end if
     call getWindow(E0, num, tmpEig, M, oldEig) 

     num = cnt_high
     HMAT(1:num,1:num) = HMAT_old(1:num,1:num) 
     call DSYEV('N', 'U', num, HMAT, M_MAX, tmpEIG, WORK, LWORK, PISTF_Conv)

     if (PISTF_CONV /= 0) then       ! Lapack error   
         PISTF_CONV = - num;       return                
     end if
     call getWindow(E0, num, tmpEig, M, Eig) 

     select case (nType)
     case (:0)
          RES = getDiffMax(M,eig, oldEig)        

     case (1:)
          RES = getSepMax(M, oldEig, cnt_high, tmpEig)
     end select

     print *, 'Cnt1=',CNT_HIGH, ' Cnt2=',(cnt_High-stepEig),' RES=', RES

     if ( RES <= ETOL ) then
          PISTF_CONV = cnt_high;       exit
     end if    

     PISTF_CONV = - cnt_high 

     cnt_low  = cnt_high + 1
     if (cnt_low >= M_MAX) then
           pistf_conv = - M_MAX;     exit
      end if

      cnt_high = cnt_high + stepCNT

      if (cnt_high > M_MAX)  cnt_high = M_MAX

  end do

  if (PISTF_CONV>0) then  ! store data
     open(99, FILE=filename, status='Replace', FORM='UNFORMATTED')
     write(99) PISTF_CONV, N
     write(99) tmpEig(1:PISTF_CONV), HMat_old(1:PISTF_CONV,1:PISTF_CONV)
     write(99) VJ(1:N, 1:PISTF_CONV)
     close(99)
  end if

  call Reorder('A', M, eig)
  
end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
