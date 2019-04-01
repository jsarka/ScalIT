!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A fast version but difficult to implement                  c
!c         It will reuse the previous generated VJ, and Hij               c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function Lan_Conv(E0, ETOL, nType, N, X, startCnt, stepCnt, stepEig,    &
                          M, M_MAX, HX, EIG, RES)
  implicit none
  double precision, intent(in)  :: E0,ETOL
  integer, intent(IN)           :: nType, N, M, M_MAX
  integer, intent(IN)           :: startCnt, stepCnt, stepEig
  double precision, intent(in)  :: X(N)  
  double precision, intent(out) :: EIG(M)   
  external                      :: HX    
  double precision, intent(out) :: RES             

!cccccccccccccccccccccccccccccccc
  double precision :: VJ(N, M_MAX+1)  
  double precision :: Alpha(M_MAX), ALPHA_OLD(M_MAX)
  double precision :: BETA(M_MAX), BETA_OLD(M_MAX)
  double precision :: oldEig(M)  
  
  integer   :: cnt_high, cnt_low, cnt_start
  integer   :: I, J, NUM, INFO

!cccccccccccccccccccccccccccccccc  
  double precision  :: getDiffMax, getSepMax
  integer           :: LANH0
  double precision  :: tmp

!ccccccccccccc   Parameters to call SSTEV subroutine cccccccccccccc
  integer, parameter  :: LDZ = 1   
  double precision    :: WORK(LDZ), Z(LDZ)

!************************************ 
  LAN_CONV = 0
  if ((stepEig >= M_MAX) .or. (startCnt > M_MAX) .or. (M > M_MAX)) &
      return

  ALPHA(1:M_MAX) = 0.0D0;   BETA(1:M_MAX)  = 0.0D0

  cnt_low   = 1;     Cnt_high  = max(startCnt, M, stepEig+1)   

  do  while ( cnt_high <= M_MAX )     

     NUM = LANH0(N, X, VJ, CNT_LOW, CNT_HIGH, HX, ALPHA, BETA)
     if (NUM <= 0) then
         LAN_CONV = -(CNT_LOW-1);      return
     end if

     ALPHA_OLD(1:CNT_HIGH) = ALPHA(1:CNT_HIGH)
     BETA_OLD (1:CNT_HIGH) = BETA (1:CNT_HIGH)
    
     num = cnt_high - stepEig
     call DSTEV('N', num, ALPHA, BETA, Z, LDZ, WORK, INFO)
     if (INFO /= 0) then       ! Lapack error   
         LAN_CONV = -num;      return                
     end if
     call getWindow(E0, num, alpha,  M, oldEig )

     num = cnt_high
     ALPHA(1:num) = ALPHA_OLD(1:num)
     BETA (1:num) = BETA_OLD (1:num)
     call DSTEV('N', num, ALPHA, BETA, Z, LDZ, WORK, INFO)
     if (INFO /= 0) then       ! Lapack error   
         LAN_CONV = -num;      return                
     end if
     call getWindow(E0, num, alpha, M, eig)

     select case (nType)
     case (:0)
          RES = getDiffMax(M,eig, oldEig)        

     case (1:)
          RES = getSepMax(M, oldEig, cnt_high, ALPHA)
     end select

     print *, 'CNT1=',CNT_HIGH, ' Cnt2=',(cnt_High-stepEig),' RES=', RES
     print *, 'Eig:', eig(1:M)

     if ( RES < ETOL ) then
          LAN_CONV = cnt_high;     exit
     end if    

     ALPHA(1:CNT_HIGH) = ALPHA_OLD(1:CNT_HIGH)
     BETA (1:CNT_HIGH) = BETA_OLD (1:CNT_HIGH)

     LAN_CONV = - cnt_high 

     cnt_low  = cnt_high + 1

     if (cnt_low >= M_MAX) then
           lan_conv = - M_MAX;      exit
      end if

      cnt_high = cnt_high + stepCnt

      if (cnt_high > M_MAX)  cnt_high = M_MAX

  end do

  call Reorder('A', M, eig)
  
end function LAN_CONV

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A Slow version but easy to implement                       c
!c         It will reuse the previous generated VJ, and Hij               c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function LAN_CONVERG(E0, ETOL, nType, N, X, startCnt, stepCnt, stepEig, &
                          M, M_MAX, HX, EIG, RES)
  implicit none
  double precision, intent(in) :: E0,ETOL
  integer, intent(IN)          :: nType, N, M, M_MAX
  integer, intent(IN)          :: startCnt, stepCnt, stepEig
  double precision, intent(in)  :: X(N)  
  double precision, intent(out) :: EIG(M)   
  external                      :: HX    
  double precision, intent(out) :: RES  

!cccccccccccccccccccccccccccccccc
  double precision :: tmp1Eig(M_MAX),tmp2Eig(M_MAX)
  double precision :: oldEig(M)  
  integer          :: cnt1, cnt2, num

!cccccccccccccccccccccccccccccccc
  integer           :: LANEIG
  double precision  :: getDiffMax, getSepMax
   
!************************************   
  LAN_CONVERG = 0
  if ((stepEig >= M_MAX) .or. (startCnt > M_MAX) .or. (M > M_MAX)) &
      return

  Cnt1  = max(startCnt, M, stepEig+1) 

  LAN_CONVERG = - M_MAX
  do  while ( CNT1 <= M_MAX )

      CNT2    = CNT1 - stepEig

      num = LANEIG(N,X,CNT1,CNT2,HX,tmp1EIG,tmp2Eig)  

      if (num <= 0) then
          LAN_CONVERG = -CNT1;         return
      end if
      
     call getWindow(E0, CNT2, tmp2Eig,  M, oldEig )
     call getWindow(E0, cnt1, tmp1Eig,  M, eig)

     select case (nType)
     case (:0)
          RES = getDiffMax(M,eig, oldEig)        

     case (1:)
           RES = getSepMax(M, oldEig, cnt1, tmp1Eig)
     end select      

     print *, 'CNT1=', CNT1,' CNT2=', CNT2, ' RES=', RES

     if ( RES < ETOL ) then
          LAN_CONVERG = CNT1;     exit
     end if    

     CNT1 = CNT1 + stepCnt

  end do

  call Reorder('A', M, eig)    
  
end 


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    LANCZOS algorithm  to get the eigen values near energy E     c
!c       Apply for Symmetric matrix H                              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function LAN(N, X, M, HX, EIG)
  implicit none
  integer, intent(IN)  :: N, M
  double precision, intent(in) :: X(N)  
  external             :: HX

  double precision, intent(out) :: EIG(M)

!cccccccccccccccccccc Local variables
  double precision :: BETA(M)  
  integer :: LANHIJ
  integer, parameter :: LDZ = 1  
  double precision   :: WORK(LDZ), Z(LDZ)    
 
  LAN = LANHIJ(N, X, M, HX, EIG, BETA)

  if ( LAN > 0 ) then
      call DSTEV('N', M, EIG, BETA, Z, LDZ, WORK,LAN)

     if (LAN == 0) then      
         LAN = M
     else
         LAN = -M     ! Lapack error
     end if
  end if

end 

!**************************************************************************

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    LANCZOS algorithm  to get the eigen values near energy E     c
!c       For two sets of, Apply for Symmetric matrix H             c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function LANEIG(N, X, M1, M2, HX, EIG1, EIG2)
  implicit none
  integer, intent(IN)  :: N, M1, M2
  double precision, intent(in) :: X(N) 
  external             :: HX
  double precision, intent(out) :: EIG1(M1)
  double precision, intent(out) :: EIG2(M2)

!ccccccccccccccc  Local Variables cccccccccccccccccccccccccc  
  integer  :: M_MAX, INFO
  double precision, allocatable :: ALPHA(:), BETA(:),TMPB(:)
  integer, parameter :: LDZ = 1
  double precision   :: WORK(LDZ), Z(LDZ)
  integer :: LANHIJ

  M_MAX = max(M1, M2);  LANEIG = 0

  allocate( ALPHA(M_MAX),BETA(M_MAX), TMPB(M_MAX), STAT=INFO)
  if (INFO /= 0)    return 

  LANEIG = LANHIJ(N, X, M_MAX, HX, ALPHA, BETA)

  if ( LANEIG > 0 ) then
      EIG1(1:M1) = ALPHA(1:M1)
      TMPB(1:M1) = BETA (1:M1)
      call DSTEV('N', M1, EIG1, TMPB, Z, LDZ, WORK, INFO)      

      if (INFO /= 0) then
           LANEIG = - M_MAX
      else
           EIG2(1:M2) = ALPHA(1:M2)
           TMPB(1:M2) = BETA(1:M2)
           call DSTEV('N',M2, EIG2, TMPB, Z, LDZ, WORK, INFO)             
           if (INFO /= 0)   LANEIG = - M_MAX
      end if
  end if

  deallocate(ALPHA, BETA, TMPB, STAT=INFO)

end function LANEIG

!**************************************************************************

