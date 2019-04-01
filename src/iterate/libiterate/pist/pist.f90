!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A fast version but difficult to implement                  c
!c         It will reuse the previous generated VJ, and Hij               c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function PIST_CONV(E0, ETOL, NTYPE, N, X, STARTCNT,  &
          STEPCNT, STEPEIG, M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  double precision, intent(in)  :: E0, ETOL
  integer, intent(IN)           :: nType, N
  double precision, intent(in)  :: X(N) 
  integer, intent(IN)           :: startCNT, stepCNT, stepEIG
  integer, intent(IN)           :: M, M_MAX
  external                      :: H0X
  integer, external             :: LinSolv   
  double precision, intent(out) :: EIG(M), RES   

!cccccccccccccccccccccccccccccccc
  double precision :: VJ(N, M_MAX)
  double precision :: HMAT(M_MAX, M_MAX), HMAT_OLD(M_MAX, M_MAX)
  double precision :: tmpEig(M_MAX)
  double precision :: oldEig(M)  
  
  integer   :: cnt_high, cnt_low, cnt_start
  integer   :: I, J, NUM

!cccccccccccccccccccccccccccccccc  
  double precision  :: getDiffMax, getSepMax
  integer           :: PISTH0

!ccccccccccccc   Parameters to call DSYEV subroutine cccccccccccccc
  integer          :: LWORK    
  double precision :: WORK (3*M_MAX) 
  LWORK = 3*M_MAX  

!************************************ 
  if ((stepEig >= M_MAX) .OR. (startCnt > M_MAX) .OR. (M > M_MAX)) then
      PIST_CONV = 0;   return
  end if

  cnt_low = 1;    cnt_high = max(startCnt, M, stepEig+1)   

  do  while ( cnt_high <= M_MAX )

     PIST_CONV = PISTH0(N,X,M_MAX,CNT_LOW,CNT_HIGH,H0X,LinSolv,VJ,HMAT_OLD)

     if (PIST_CONV < 0) then
         PIST_CONV = - CNT_HIGH;  return
     end if    

     num = cnt_high-stepEig
     HMAT(1:num,1:num) = HMAT_OLD(1:num,1:num)     
     call DSYEV('N', 'U', num, HMAT, M_MAX, tmpEIG, WORK, LWORK, PIST_Conv)
     if (PIST_CONV /= 0) then       ! Lapack error   
         PIST_CONV = -num;    return                
     end if
     call getWindow(E0, num, tmpEig, M, oldEig) 

     num = cnt_high
     HMAT(1:num,1:num) = HMAT_old(1:num,1:num) 
     call DSYEV('N', 'U', num, HMAT, M_MAX, tmpEIG, WORK, LWORK, PIST_Conv)

     if (PIST_CONV /= 0) then       ! Lapack error   
         PIST_CONV = - num;     return                
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
          PIST_CONV = cnt_high;    exit
     end if    

     PIST_CONV = - cnt_high 

     cnt_low  = cnt_high + 1
     if (cnt_low >= M_MAX) then
         pist_conv = - M_MAX;   exit
      end if

      cnt_high = cnt_high + stepCNT

      if (cnt_high > M_MAX) cnt_high = M_MAX

  end do

  call Reorder('A', M, eig)
  
end function PIST_CONV


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A slow version but easy to implement                       c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function PIST_CONVERG(E0, ETOL, NTYPE, N, X, STARTCNT, &
             STEPCNT, STEPEIG, M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  double precision, intent(in)  :: E0, ETOL
  integer, intent(IN        )   :: nType, N
  double precision,intent(in)   :: X (N) 
  integer, intent(IN)           :: startCNT, stepCNT, stepEIG
  integer, intent(IN)           :: M, M_MAX
  external                      :: H0X
  integer, external             :: LinSolv
  double precision, intent(out) :: EIG(M)
  double precision, intent(out) :: RES 

!cccccccccccccccccccccccccccccccc
  double precision :: tmp1Eig(M_MAX)
  double precision :: tmp2Eig(M_MAX)
  double precision :: oldEig(M) 
  integer          :: cnt1, cnt2

!cccccccccccccccccccccccccccccccc
  integer          :: PIST_EIG
  double precision :: getDiffMax, getSepMax
   
!************************************   

  if ((stepEig >= M_MAX) .OR. (startCnt > M_MAX) .OR. (M > M_MAX)) then
      PIST_CONVERG = 0;     return
  end if

  Cnt1  = max(startCnt, M, stepEig+1) 

  PIST_CONVERG = - M_MAX
  do  while ( CNT1 <=M_MAX )
      CNT2    = CNT1 - stepEig

      PIST_CONVERG = PIST_EIG(N,X,CNT1,CNT2,H0X,LinSolv,tmp1EIG,tmp2Eig)

      if (PIST_CONVERG <= 0) then
          PIST_CONVERG = -CNT1;   return
      end if
       
      call getWindow(E0, CNT1, tmp1Eig, M, Eig) 
      call getWindow(E0, CNT2, tmp2Eig, M, oldEig)     
      
      select case (nType)
      case (:0)
          RES = getDiffMax(M, eig, oldEig)        

      case (1:)
           RES = getSepMax(M, oldEig, cnt1, tmp1Eig)
      end select  

      print *, 'Cnt1=', CNT1,' Cnt2=', CNT2, ' RES=', RES
  
      if ( RES < ETOL ) then
          PIST_CONVERG = CNT1;   exit
      end if    

      CNT1 = CNT1 + stepCnt

      if (CNT1>=M_MAX) CNT1=M_MAX

  end do

  call Reorder('A', M, eig)
  
end function PIST_CONVERG


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    PIST algorithm  to get the eigen values near energy E        c
!c             Apply for Symmetric matrix H                        c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function PIST(N, X, M, H0X, LinSolv, EIG)
  implicit none
  integer, intent(IN)  :: N, M
  double precision,intent(in) :: X(N)  
  external             :: H0X
  integer, external    :: LinSolv
  double precision,intent(out) :: EIG(M)

!cccccccccccccccccccc Local variables
  integer  :: PISTHIJ
  integer  :: LWORK  
  double precision :: WORK (3*M) 
  double precision :: HMAT(M,M)  ! HMAT = <U|H|U>

  LWORK = 3*M     

  PIST = PISTHIJ(N, X, M, H0X, LinSolv, HMAT)

  if ( PIST > 0 ) then
      call DSYEV('N', 'U', M, HMAT, M, EIG, WORK, LWORK, PIST)

     if (PIST == 0) then      
         PIST = M
     else
         PIST = -M    
     end if
  end if

end function PIST
!********************************************************************


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    PIST algorithm  to get the eigen values near energy E        c
!c     For two sets of Eigs.  Apply for Symmetric matrix H         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function PIST_EIG(N, X,  M1, M2, H0X, LinSolv, EIG1, EIG2)
  implicit none
  integer, intent(IN)  :: N                       
  double precision, intent(in) :: X(N)
  integer, intent(IN)  :: M1, M2                  
  external             :: H0X
  integer, external    :: LinSolv

  double precision, intent(out) :: EIG1(M1)
  double precision, intent(out) :: EIG2(M2)

!ccccccccccccccc  Local Variables cccccccccccccccccccccccccc  
  integer  :: M_MAX, INFO, LWORK
  double precision, allocatable :: HMAT(:,:), HOLD(:,:)  
  double precision, allocatable :: WORK(:)
  integer :: PISTHIJ

  M_MAX = max(M1, M2);   LWORK = 3*M_MAX  
  PIST_EIG = -1

  allocate( HMAT(M_MAX, M_MAX),HOLD(M_MAX, M_MAX),        &
            WORK(LWORK), STAT=INFO)
  if (INFO /= 0) return 

  PIST_EIG = PISTHIJ(N, X, M_MAX, H0X, LinSolv, HMAT)

  if ( PIST_EIG > 0 ) then
      HOLD(1:M_MAX, 1:M_MAX) = HMAT(1:M_MAX, 1:M_MAX)
      call DSYEV('N', 'U', M1, HMAT, M_MAX, EIG1, WORK, LWORK, INFO)

      if (INFO /= 0) then
           PIST_EIG = - M_MAX
      else
           call DSYEV('N', 'U', M2, HOLD, M_MAX, EIG2, WORK, LWORK, INFO)             
           if (INFO /= 0)  PIST_EIG = - M_MAX           
      end if
  end if

  deallocate( HOLD, HMAT, WORK, STAT=INFO)

end function PIST_EIG
!**************************************************************************
