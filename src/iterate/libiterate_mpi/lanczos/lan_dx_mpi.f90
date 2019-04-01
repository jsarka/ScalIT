!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A fast version but difficult to implement                  c
!c         It will reuse the previous generated VJ, and Hij               c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function lan_Conv_DX_MPI( MYID, ROOTID, E0, ETOL, nType, N, X,    &
                       startCnt, stepCnt, stepEig, M, M_MAX, HX, EIG, RES)
  implicit none
  include 'mpif.h'
  integer, intent(IN)          :: MYID, ROOTID
  double precision, intent(in) :: E0,ETOL
  integer, intent(IN)          :: nType, N, M, M_MAX
  integer, intent(IN)          :: startCnt, stepCnt, stepEig
  double complex,  intent(in)  :: X(N) 
  double precision,intent(out) :: EIG(M)
  external                     :: HX    
  double precision, intent(out) :: RES             


!cccccccccccccccccccccccccccccccc
  double complex   :: VJ(N, M_MAX+1)     
  double precision :: ALPHA(M_MAX), ALPHA_OLD(M_MAX)
  double precision :: BETA(M_MAX), BETA_OLD(M_MAX)
  double precision :: oldEig(M)
  
  integer   :: cnt_high, cnt_low, cnt_start
  integer   :: I, J, NUM, ierr,INFO(2)

!cccccccccccccccccccccccccccccccc  
  double precision  :: getDiffMax, getSepMax
  integer           :: LANH0_dx_MPI
  double precision  :: tmp

!ccccccccccccc   Parameters to call SSTEV subroutine cccccccccccccc
  integer, parameter   :: LDZ = 1   
  double precision, dimension(LDZ) :: WORK, Z

!************************************ 
  LAN_CONV_dx_MPI = 0
  if ((stepEig >= M_MAX) .or. (startCnt > M_MAX) .or. (M > M_MAX))   return

  ALPHA(1:M_MAX) = 0.0D0
  BETA(1:M_MAX)  = 0.0D0

  cnt_low   = 1
  Cnt_high  = max(startCnt, M, stepEig+1)   

  LAN_CONV_dx_MPI = -Cnt_high

  do  while ( cnt_high <= M_MAX )     

     NUM = LANH0_dx_MPI( MYID, ROOTID, N, X, VJ, CNT_LOW,   &
                        CNT_HIGH, HX, ALPHA, BETA)
     if (NUM <= 0) then
         LAN_CONV_dx_MPI = -(CNT_LOW-1)
         return
     end if

     ALPHA_OLD(1:CNT_HIGH) = ALPHA(1:CNT_HIGH)
     BETA_OLD (1:CNT_HIGH) = BETA (1:CNT_HIGH)

     if (myID == rootID) then   

        num = cnt_high - stepEig

        call DSTEV('N', num, ALPHA, BETA, Z, LDZ, WORK, INFO(1))

        if (INFO(1) == 0) then       
            call getWindow(E0, num, alpha, M, oldeig)

            num = cnt_high 
            ALPHA(1:num) = ALPHA_OLD(1:num)
            BETA (1:num) = BETA_OLD (1:num)
            call DSTEV('N', num, ALPHA, BETA, Z, LDZ, WORK, INFO(2))

            if (INFO(2) == 0) then       
                call getWindow(E0, num, alpha,  M, eig )

                select case (nType)
                case (:0)
                     RES = getDiffMax(M,eig, oldEig)        

                case (1:)
                     RES = getSepMax(M, oldEig, cnt_high, alpha)
                end select

                print *, 'CNT1=',CNT_HIGH, ' Cnt2=',(CNT_HIGH-STEPEIG),' RES=', RES
                print *, 'Eig:', eig
            end if
        end if

     end if

     ALPHA(1:CNT_HIGH) = ALPHA_OLD(1:CNT_HIGH)
     BETA (1:CNT_HIGH) = BETA_OLD (1:CNT_HIGH)

     call MPI_BCAST(INFO, 2, MPI_INTEGER, ROOTID, MPI_COMM_WORLD, IERR)
      
     if ((INFO(1)==0) .and. (INFO(2)==0)) then
         call MPI_BCAST(RES, 1, MPI_DOUBLE_PRECISION, ROOTID, MPI_COMM_WORLD, IERR)
     else
         LAN_CONV_dx_MPI = -num
         return 
     end if

     if ( RES < ETOL ) then
          LAN_CONV_dx_MPI = cnt_high
          exit
     end if    

     LAN_CONV_dx_MPI = - cnt_high 

     cnt_low  = cnt_high + 1
     if (cnt_low >= M_MAX) then
           lan_conv_dx_mpi = - M_MAX
           exit
      end if
      cnt_high = cnt_high + stepCnt
      if (cnt_high > M_MAX) then
           cnt_high = M_MAX
      end if
  end do

  if (myID==rootID)  call Reorder('A', M, eig)
  
end 

!****************************************************************************
double precision function CRP_lan_Conv_DX_MPI( MYID, ROOTID, ETOL, N, X,    &
                       startCnt, stepCnt, M_MAX, HX, RES)  ! EIG, RES
  implicit none
  include 'mpif.h'
  integer, intent(IN)          :: MYID, ROOTID
  double precision, intent(in) :: ETOL
  integer, intent(IN)          :: N, M_MAX
  integer, intent(IN)          :: startCnt, stepCnt
  double complex, intent(in)   :: X(N)
  external                     :: HX    
  double precision, intent(out):: RES             


!cccccccccccccccccccccccccccccccc
  double complex   :: VJ(N, M_MAX+1)    
  double precision :: ALPHA(M_MAX), ALPHA_OLD(M_MAX)
  double precision :: BETA(M_MAX),  BETA_OLD(M_MAX)
  
  integer   :: cnt_high, cnt_low, cnt_start
  integer   :: NUM, ierr, INFO
  double precision   ::  dINFO(3)

!cccccccccccccccccccccccccccccccc  
  integer           :: LANH0_dx_MPI
  double precision  :: tmp

!ccccccccccccc   Parameters to call SSTEV subroutine cccccccccccccc
  integer, parameter   :: LDZ = 1   
  double precision, dimension(LDZ) :: WORK, Z

!************************************ 
  CRP_LAN_CONV_dx_MPI = 0.0D0
  if (startCnt > M_MAX)    return

  ALPHA(1:M_MAX) = 0.0D0
  BETA(1:M_MAX)  = 0.0D0

  cnt_low   = 1
  Cnt_high  = startCnt   

  do  while ( cnt_high <= M_MAX )     

     NUM = LANH0_dx_MPI( MYID, ROOTID, N, X, VJ, CNT_LOW,   &
                        CNT_HIGH, HX, ALPHA, BETA)
     if (NUM <= 0) return

     ALPHA_OLD(1:CNT_HIGH) = ALPHA(1:CNT_HIGH)
     BETA_OLD (1:CNT_HIGH) = BETA (1:CNT_HIGH)

     if (myID == rootID) then   

        num = cnt_high 

        call DSTEV('N', num, ALPHA, BETA, Z, LDZ, WORK, INFO)

        if (INFO == 0) then       
           RES = sum(alpha(1:num)) - CRP_LAN_CONV_DX_MPI
           CRP_LAN_CONV_DX_MPI = CRP_LAN_CONV_DX_MPI + RES
           RES = DABS(RES)
           print *, 'CNT1=',CNT_HIGH, ' CRP=', CRP_LAN_CONV_DX_MPI,' Error=', RES
           dInfo(3) = 1.0D0
        else
           dInfo(3) = 0.0D0
        end if 
        dInfo(1) = RES
        dInfo(2) = CRP_LAN_CONV_DX_MPI
     end if

     ALPHA(1:CNT_HIGH) = ALPHA_OLD(1:CNT_HIGH)
     BETA (1:CNT_HIGH) = BETA_OLD (1:CNT_HIGH)

     call MPI_BCAST(dINFO, 3, MPI_DOUBLE_PRECISION, ROOTID, MPI_COMM_WORLD, IERR)
      
     RES = dInfo(1)
     CRP_LAN_CONV_DX_MPI = dInfo(2)

     if (DABS(dINFO(1)) < 1.0D-5) return

     if ( RES < ETOL )  return

     cnt_low  = cnt_high + 1
     if (cnt_low >= M_MAX) return

      cnt_high = cnt_high + stepCnt
      if (cnt_high > M_MAX)    cnt_high = M_MAX
    
  end do
  
end 


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    LANCZOS algorithm  to get the eigen values near energy E     c
!c       Apply for Symmetric matrix H                              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Get the eigen values via Full Reorthogonization Lanczos:
!
! integer function Lan(N, X, M, HX, EIG)
!
!
integer function LAN_dx_MPI(MYID, ROOTID, N, X, M, HX, EIG)
  implicit none
  include 'mpif.h'
  integer, intent(IN)  :: MYID, ROOTID
  integer, intent(IN)  :: N, M
  double complex, intent(in) :: X(N)
  external             :: HX
  double precision,intent(out) :: EIG(M)
 

!cccccccccccccccccccc Local variables
  double precision :: BETA(M)
  integer :: LANHIJ_dx_MPI
  integer, parameter :: LDZ = 1  
  double precision, dimension(1) :: WORK, Z    
 
  LAN_dx_MPI = LANHIJ_dx_MPI( MYID, ROOTID, N, X, M, HX, EIG, BETA)

  if ( LAN_dx_MPI > 0 ) then
     if (myid == rootID) then
          call DSTEV('N', M, EIG, BETA, Z, LDZ, WORK, LAN_dx_MPI)
     end if

     if (LAN_dx_MPI == 0) then      
         LAN_dx_MPI = M
     else
         LAN_dx_MPI = -M     ! Lapack error
     end if
  end if

  call MPI_BCAST(LAN_dx_MPI, 1, MPI_INTEGER, ROOTID, MPI_COMM_WORLD)

end function LAN_dx_MPI

!**************************************************************************


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    LANCZOS algorithm  to get the eigen values near energy E     c
!c    For two sets og eigenvalues, Apply for Symmetric matrix H    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function LANEIG_dx_MPI( MYID, ROOTID, N, X, M1, M2,  &
                               HX, EIG1, EIG2)
  implicit none
  include 'mpif.h'
  integer, intent(IN)  :: MYID, ROOTID
  integer, intent(IN)  :: N, M1, M2
  double complex, intent(in) :: X(N)
  external             :: HX
  double precision, intent(out) :: EIG1(M1)
  double precision, intent(out) :: EIG2(M2)


!ccccccccccccccc  Local Variables cccccccccccccccccccccccccc  
  integer  :: M_MAX, INFO, ierr
  double precision, allocatable :: ALPHA(:), BETA(:),TMPB(:)
  integer, parameter :: LDZ = 1
  double precision, dimension(LDZ) :: WORK, Z
  integer :: LANHIJ_dx_MPI

  M_MAX = max(M1, M2);  LANEIG_dx_MPI = 0

  allocate( ALPHA(M_MAX),BETA(M_MAX), TMPB(M_MAX), STAT=INFO)
  if (INFO /= 0)    return 

  LANEIG_dx_MPI = LANHIJ_dx_MPI( MYID, ROOTID, N, X, M_MAX,  &
                                HX, ALPHA, BETA)

  if ( LANEIG_dx_MPI > 0 ) then
    if (MYID==ROOTID) then
      EIG1(1:M1) = ALPHA(1:M1)
      TMPB(1:M1) = BETA (1:M1)
      call DSTEV('N', M1, EIG1, TMPB, Z, LDZ, WORK, INFO)      

      if (INFO /= 0) then
           LANEIG_dx_MPI = - M_MAX
      else
           EIG2(1:M2) = ALPHA(1:M2)
           TMPB(1:M2) = BETA(1:M2)
           call DSTEV('N',M2, EIG2, TMPB, Z, LDZ, WORK, INFO)             
           if (INFO /= 0)     LANEIG_dx_MPI = - M_MAX           
      end if
    end if
  end if

  call MPI_BCAST(LANEIG_dx_MPI, 1, MPI_INTEGER, ROOTID, MPI_COMM_WORLD, ierr)

  deallocate(ALPHA, BETA, TMPB, STAT=INFO)

end 

!**************************************************************************


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A Slow version but easy to implement                       c
!c         It will reuse the previous generated VJ, and Hij               c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function LAN_CONVERG_dx_MPI( MYID, ROOTID, E0, ETOL, nType,  &
                   N, X, startCnt, stepCnt, stepEig, M, M_MAX, HX, EIG, RES)
  implicit none
  include 'mpif.h'
  integer, intent(IN)          :: MYID, ROOTID
  double precision, intent(in) :: E0,ETOL
  integer, intent(IN)          :: nType, N, M, M_MAX
  integer, intent(IN)          :: startCnt, stepCnt, stepEig
  double complex,  intent(in)  :: X(N)  
  double precision, intent(out) :: EIG(M)   
  external                      :: HX    
  double precision, intent(out) :: RES  


!cccccccccccccccccccccccccccccccc
  double precision :: tmp1Eig(M_MAX),tmp2Eig(M_MAX),oldEig(M)
  integer          ::  cnt1, cnt2, num, IERR

!cccccccccccccccccccccccccccccccc
  integer           :: LANEIG_dx_MPI
  double precision  :: getDiffMax, getSepMax
   
!************************************   
  LAN_CONVERG_dx_MPI = 0
  if ((stepEig >= M_MAX) .or. (startCnt > M_MAX) .or. (M > M_MAX))  return

  Cnt1  = max(startCnt, M, stepEig+1) 

  LAN_CONVERG_dx_MPI = - CNT1
  do  while ( CNT1 <=M_MAX )

      CNT2    = CNT1 - stepEig

      num = LANEIG_dx_MPI( MYID,ROOTID,N,X,CNT1,CNT2,HX,tmp1EIG,tmp2Eig)  

      if ( num <= 0) then
          LAN_CONVERG_dx_MPI = - CNT1
          return
      end if
      
     if (myID == rootID) then
          call getWindow(E0, CNT2, tmp2Eig,  M, oldEig )
          call getWindow(E0, cnt1, tmp1Eig,  M, eig)

          select case (nType)
          case (:0)
               RES = getDiffMax(M,eig, oldEig)        

          case (1:)
               RES = getSepMax(M, oldEig, cnt1, tmp1Eig)
          end select      

          print *, 'CNT1=', CNT1,' CNT2=', CNT2, ' RES=', RES
     end if
 
     call MPI_BCAST(RES, 1, MPI_DOUBLE_PRECISION, rootID, MPI_COMM_WORLD, IERR)

     if ( RES < ETOL ) then
          LAN_CONVERG_dx_MPI = CNT1
          exit
     end if    

     CNT1 = CNT1 + stepCnt
  end do

  if (myID == rootID)       call Reorder('A', M, eig)  
 
end function LAN_CONVERG_dx_MPI

!**************************************************************************
double precision function CRP_LAN_CONVERG_dx_MPI( MYID, ROOTID, ETOL,    &
                   N, X, startCnt, stepCnt, M_MAX, HX, RES)
  implicit none
  include 'mpif.h'
  integer, intent(IN)          :: MYID, ROOTID
  double precision, intent(in) :: ETOL
  integer, intent(IN)          :: N, M_MAX
  integer, intent(IN)          :: startCnt, stepCnt
  double complex, intent(in)   :: X(N)
  external                     :: HX    
  double precision, intent(out):: RES  


!cccccccccccccccccccccccccccccccc
  integer          :: cnt1, cnt2, num, IERR
  double precision :: EIG(M_MAX)
  integer           :: LAN_dx_MPI
   
!************************************   
  CRP_LAN_CONVERG_dx_MPI = 0.0D0
  if (startCnt > M_MAX)  return

  Cnt1  = startCnt

  do  while ( CNT1 <=M_MAX )

      num = LAN_dx_MPI( MYID, ROOTID, N, X, CNT1, HX, EIG)  

      if ( num <= 0)  return
      
      if (myID == rootID) then
          RES = sum(EIG(1:CNT1))-CRP_LAN_CONVERG_DX_MPI
          CRP_LAN_CONVERG_DX_MPI = RES + CRP_LAN_CONVERG_DX_MPI
          RES = DABS(RES)

          print *, 'CNT1=', CNT1, ' RES=', RES, ' CRP=',CRP_LAN_CONVERG_DX_MPI
     end if
 
     call MPI_BCAST(RES, 1, MPI_DOUBLE_PRECISION, rootID, MPI_COMM_WORLD,IERR)

     if ( RES < ETOL ) return

     CNT1 = CNT1 + stepCnt
  end do

end function

!**************************************************************************
