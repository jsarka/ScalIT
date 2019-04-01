!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A fast version but difficult to implement                  c
!c         It will reuse the previous generated VJ, and Hij               c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


integer function LAN_CONV_DX(E0, ETOL, nType, N, X, startCnt, stepCnt,    &
                         stepEig, M, M_MAX, HXCX, EIG, RES)
  implicit none
  double precision, intent(in) :: E0,ETOL
  integer, intent(IN)          :: nType, N, M, M_MAX
  integer, intent(IN)          :: startCnt, stepCnt, stepEig
  double complex,  intent(in)  :: X(N)  
  double precision,intent(out) :: EIG(M)   
  external                     :: HXCX    
  double precision, intent(out):: RES             

!cccccccccccccccccccccccccccccccc
  double complex :: VJ(N, M_MAX+1)
  double precision :: ALPHA(M_MAX), ALPHA_OLD(M_MAX)
  double precision :: BETA(M_MAX), BETA_OLD(M_MAX)
  double precision :: oldEig(M)  
  
  integer :: cnt_high, cnt_low, cnt_start
  integer :: I, J, NUM, INFO

!cccccccccccccccccccccccccccccccc  
  double precision :: getDiffMax, getSepMax
  integer          :: LANH0_DX
  double precision :: tmp

!ccccccccccccc   Parameters to call SSTEV subroutine cccccccccccccc
  integer, parameter :: LDZ = 1   
  double precision   :: WORK(LDZ), Z(LDZ)

!************************************ 
  LAN_CONV_DX = 0
  if ((stepEig >= M_MAX) .or. (startCnt > M_MAX) .or. (M > M_MAX)) &
      return

  ALPHA(1:M_MAX) = 0.0D0;  BETA(1:M_MAX)  = 0.0D0

  cnt_low   = 1
  Cnt_high  = max(startCnt, M, stepEig+1)   

  do  while ( cnt_high <= M_MAX )     

     NUM = LANH0_DX(N, X, VJ, CNT_LOW, CNT_HIGH, HXcx, ALPHA, BETA)
     if (NUM <= 0) then
         LAN_CONV_DX = -(CNT_LOW-1);         return
     end if

     ALPHA_OLD(1:CNT_HIGH) = ALPHA(1:CNT_HIGH)
     BETA_OLD(1:CNT_HIGH) = BETA(1:CNT_HIGH)
    
     num = cnt_high - stepEig
     call DSTEV('N', num, ALPHA, BETA, Z, LDZ, WORK, INFO)
     if (INFO /= 0) then       ! Lapack error   
         LAN_CONV_DX = -num;          return                
     end if
     call getWindow(E0, num, alpha, M, oldEig)

     num = cnt_high
     ALPHA(1:num) = ALPHA_OLD(1:num)
     BETA(1:num) = BETA_OLD(1:num)
     call DSTEV('N', num, ALPHA, BETA, Z, LDZ, WORK, INFO)
     if (INFO /= 0) then       ! Lapack error   
         LAN_CONV_DX = -num;       return                
     end if
     call getWindow(E0, cnt_high, alpha, M, eig)


     select case (nType)
     case (:0)
          RES = getDiffMax(M,eig, oldEig)        

     case (1:)
          RES = getSepMax(M, oldEig, cnt_high, ALPHA)
     end select

     print *, 'CNT1=',CNT_HIGH, ' Cnt2=',(cnt_High-stepEig),' RES=', RES
     print *, 'Eig:', eig(1:M)

     if ( RES < ETOL ) then
          LAN_CONV_DX = cnt_high;      exit
     end if    

     ALPHA(1:CNT_HIGH) = ALPHA_OLD(1:CNT_HIGH)
     BETA(1:CNT_HIGH) = BETA_OLD(1:CNT_HIGH)

     LAN_CONV_DX = - cnt_high 

     cnt_low  = cnt_high + 1
     if (cnt_low >= M_MAX) then
           lan_conv_dx = - M_MAX;           exit
      end if
      cnt_high = cnt_high + stepCnt
      if (cnt_high > M_MAX)    cnt_high = M_MAX
  end do

  call Reorder('A', M, eig)
  
end function LAN_CONV_DX

!**********************************************************
!c    Modified Lan_Conv_DX for CRP calculation            *
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function CRP_LAN_CONV_DX(ETOL, N, X, startCnt, stepCnt, &
                         M_MAX, HXCX, RES)
  implicit none
  double precision, intent(in) :: ETOL
  integer, intent(IN)          :: N, M_MAX
  integer, intent(IN)          :: startCnt, stepCnt
  double complex,  intent(in)  :: X(N)  
  external                      :: HXCX    
  double precision, intent(out) :: RES             

!cccccccccccccccccccccccccccccccc
  double complex   :: VJ(N, M_MAX+1) 
  double precision :: ALPHA(M_MAX), ALPHA_OLD(M_MAX)
  double precision :: BETA(M_MAX), BETA_OLD(M_MAX)
  
  integer :: cnt_high, cnt_low, cnt_start
  integer ::  NUM, INFO

!cccccccccccccccccccccccccccccccc  
  integer          :: LANH0_DX
  double precision :: tmp

!ccccccccccccc   Parameters to call SSTEV subroutine cccccccccccccc
  integer, parameter :: LDZ = 1   
  double precision   :: WORK(LDZ), Z(LDZ)

!************************************ 
  CRP_LAN_CONV_DX = 0.0D0
  if (startCnt > M_MAX) return

  ALPHA(1:M_MAX) = 0.0D0;   BETA(1:M_MAX)  = 0.0D0

  cnt_low   = 1;    Cnt_high  = startCnt

  do  while ( cnt_high <= M_MAX )     

     NUM = LANH0_DX(N, X, VJ, CNT_LOW, CNT_HIGH, HXcx, ALPHA, BETA)
     if (NUM <= 0)     return

     ALPHA_old(1:CNT_HIGH) = ALPHA(1:CNT_HIGH)
     BETA_old(1:CNT_HIGH) = BETA(1:CNT_HIGH)

     num = cnt_high
     call DSTEV('N', num, ALPHA, BETA, Z, LDZ, WORK, INFO)
     if (INFO /= 0)    return                

     RES = sum(alpha(1:num)) - CRP_LAN_CONV_DX

     CRP_LAN_CONV_DX = RES + CRP_LAN_CONV_DX

     RES = DABS(RES)

     print *, 'CNT1=',CNT_HIGH, ' CRP=', CRP_LAN_CONV_DX, ' Error:', RES

     if ( RES < ETOL )   exit

     ALPHA(1:CNT_HIGH) = ALPHA_OLD(1:CNT_HIGH)
     BETA(1:CNT_HIGH) = BETA_OLD(1:CNT_HIGH)

     cnt_low  = cnt_high + 1
     if (cnt_low >= M_MAX)   exit

      cnt_high = cnt_high + stepCnt

      if (cnt_high > M_MAX)   cnt_high = M_MAX     
  end do
  
end function CRP_LAN_CONV_DX


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A Slow version but easy to implement                       c
!c         It will reuse the previous generated VJ, and Hij               c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function LAN_CONVERG_DX(E0, ETOL, nType, N, X, startCnt, stepCnt,   &
                      stepEig, M, M_MAX, HXCX, EIG, RES)
  implicit none
  double precision, intent(in) :: E0,ETOL
  integer, intent(IN)          :: nType, N, M, M_MAX
  integer, intent(IN)          :: startCnt, stepCnt, stepEig
  double complex, intent(in)   :: X(N)  
  double precision,intent(out) :: EIG(M)   
  external                     :: HXCX    
  double precision, intent(out) :: RES  

!cccccccccccccccccccccccccccccccc
  double precision :: tmp1Eig(M_MAX), tmp2Eig(M_MAX)
  double precision :: oldEig(M)  
  integer          ::  cnt1, cnt2, num

!cccccccccccccccccccccccccccccccc
  integer          :: LANEIG_DX
  double precision :: getDiffMax, getSepMax
   
!************************************   
  LAN_CONVERG_DX = 0
  if ((stepEig >= M_MAX) .or. (startCnt > M_MAX) .or. (M > M_MAX)) &
      return

  Cnt1  = max(startCnt, M, stepEig+1) 

  LAN_CONVERG_DX = - M_MAX
  do  while ( CNT1 <= M_MAX )

      CNT2    = CNT1 - stepEig

      num = LANEIG_DX(N,X,CNT1,CNT2,HXcx,tmp1EIG,tmp2Eig)  

      if (num <= 0) then
          LAN_CONVERG_DX = -CNT1;           return
      end if
      
     call getWindow(E0, CNT2, tmp2Eig, M, oldEig)
     call getWindow(E0, cnt1, tmp1Eig, M, eig)

     select case (nType)
     case (:0)
          RES = getDiffMax(M,eig, oldEig)        

     case (1:)
           RES = getSepMax(M, oldEig, cnt1, tmp1Eig)
     end select      

     print *, 'CNT1=', CNT1,' CNT2=', CNT2, ' RES=', RES

     if ( RES < ETOL ) then
          LAN_CONVERG_DX = CNT1;           exit
     end if    

     CNT1 = CNT1 + stepCnt

  end do

  call Reorder('A', M, eig)    
  
end function LAN_CONVERG_DX

!*************************************************************************
double precision function CRP_LAN_CONVERG_DX(ETOL, N, X, startCnt, stepCnt,   &
                       M_MAX, HXCX, RES)  !EIG, RES)
  implicit none
  double precision, intent(in) :: ETOL
  integer, intent(IN)          ::  N,  M_MAX
  integer, intent(IN)          :: startCnt, stepCnt
  double complex,  intent(in)  :: X(N)  
  external                     :: HXCX    
  double precision, intent(out) :: RES  

!cccccccccccccccccccccccccccccccc
  integer          ::  cnt1, cnt2, num
  double precision ::  Eig(M_MAX)

!cccccccccccccccccccccccccccccccc
  integer :: LAN_DX
   
!************************************   
  CRP_LAN_CONVERG_DX = 0

  if (startCnt > M_MAX)   return

  Cnt1  = startCnt

  CRP_LAN_CONVERG_DX = 0.0D0

  do  while ( CNT1 <= M_MAX )

      num = LAN_DX(N,X,CNT1,HXcx,EIG)  

      if (num <= 0) return

      RES = sum(EIG(1:CNT1))-CRP_LAN_CONVERG_DX
      CRP_LAN_CONVERG_DX = RES + CRP_LAN_CONVERG_DX
      RES = DABS(RES)

      print *, 'CNT1=', CNT1, ' CRP=',CRP_LAN_CONVERG_DX,'  Error=', RES

      if ( RES < ETOL ) exit

      CNT1 = CNT1 + stepCnt

  end do
  
end function CRP_LAN_CONVERG_DX

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    LANCZOS algorithm  to get the eigen values near energy E     c
!c       Apply for Symmetric matrix H                              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function LAN_DX(N, X, M, HXCX, EIG)
  implicit none
  integer, intent(IN)  :: N, M
  double complex, intent(in) :: X(N)
  external             :: HXCX

  double precision, intent(out) :: EIG(M)

!cccccccccccccccccccc Local variables
  double precision  :: BETA(M)  
  integer :: LANHIJ_DX
  integer, parameter :: LDZ = 1  
  double precision :: WORK(LDZ), Z(LDZ)    
 
  LAN_DX = LANHIJ_DX(N, X, M, HXCX, EIG, BETA)

  if ( LAN_DX > 0 ) then
      call DSTEV('N', M, EIG, BETA, Z, LDZ, WORK,LAN_DX)

     if (LAN_DX == 0) then      
         LAN_DX = M
     else
         LAN_DX = -M    
     end if
  end if

end function LAN_DX

!***********************************************************************


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    LANCZOS algorithm  to get the eigen values near energy E     c
!c       Apply for Symmetric matrix H                              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
integer function LANEIG_DX(N, X, M1, M2, HXCX, EIG1, EIG2)
  implicit none
  integer, intent(IN) :: N, M1, M2
  double complex,  intent(in) :: X(N) 
  external             :: HXCX
  double precision, intent(out) :: EIG1(M1)
  double precision, intent(out) :: EIG2(M2)

!ccccccccccccccc  Local Variables cccccccccccccccccccccccccc  
  integer :: M_MAX, INFO
  double precision, allocatable :: ALPHA(:), BETA(:),TMPB(:)
  integer, parameter :: LDZ = 1
  double precision   :: WORK(LDZ), Z(LDZ)
  integer :: LANHIJ_DX

  M_MAX = max(M1, M2);   LANEIG_DX = 0

  allocate( ALPHA(M_MAX),BETA(M_MAX), TMPB(M_MAX), STAT=INFO)
  if (INFO /= 0)    return 

  LANEIG_DX = LANHIJ_DX(N, X, M_MAX, HXCX, ALPHA, BETA)

  if ( LANEIG_DX > 0 ) then
      EIG1(1:M1) = ALPHA(1:M1)
      TMPB(1:M1) = BETA (1:M1)
      call DSTEV('N', M1, EIG1, TMPB, Z, LDZ, WORK, INFO)      

      if (INFO /= 0) then
           LANEIG_DX = - M_MAX
      else
           EIG2(1:M2) = ALPHA(1:M2)
           TMPB(1:M2) = BETA(1:M2)
           call DSTEV('N',M2, EIG2, TMPB, Z, LDZ, WORK, INFO)             
           if (INFO /= 0) then
                  LANEIG_DX = - M_MAX
           end if
      end if
  end if

  deallocate(ALPHA, BETA, TMPB, STAT=INFO)

end function LANEIG_DX

!**************************************************************************

