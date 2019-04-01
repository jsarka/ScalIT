!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A fast version but difficult to implement                  c
!c         It will reuse the previous generated VJ, and Hij               c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function lan_Conv_MPI(MYID, ROOTID, E0, ETOL, nType, N, X,    &
                       startCnt, stepCnt, stepEig, M, M_MAX, HX, EIG, RES)
  implicit none
  include 'mpif.h'
  integer, intent(IN)          :: MYID, ROOTID
  double precision, intent(in) :: E0,ETOL
  integer, intent(IN)          :: nType, N, M, M_MAX
  integer, intent(IN)          :: startCnt, stepCnt, stepEig
  double precision,intent(in)  :: X(N)
  double precision,intent(out) :: EIG(N)   
  external                     :: HX    
  double precision, intent(out):: RES             


!cccccccccccccccccccccccccccccccc
  double precision :: VJ(N, M_MAX+1)    
  double precision :: ALPHA(M_MAX), ALPHA_OLD(M_MAX)
  double precision :: BETA(M_MAX), BETA_OLD(M_MAX)
  double precision :: oldEig(M)
  
  integer   :: cnt_high, cnt_low, cnt_start
  integer   :: I, J, NUM, INFO(2)

!cccccccccccccccccccccccccccccccc  
  double precision  :: getDiffMax, getSepMax
  integer           :: LANH0_MPI, IERR
  double precision  :: tmp

!ccccccccccccc   Parameters to call SSTEV subroutine cccccccccccccc
  integer, parameter   :: LDZ = 1   
  double precision, dimension(LDZ) :: WORK, Z

!************************************ 
  LAN_CONV_MPI = 0
  if ((stepEig >= M_MAX) .or. (startCnt > M_MAX) .or. (M > M_MAX))   return

  ALPHA(1:M_MAX) = 0.0D0
  BETA(1:M_MAX)  = 0.0D0

  cnt_low   = 1
  Cnt_high  = max(startCnt, M, stepEig+1)   

  LAN_CONV_MPI = -Cnt_high

  do  while ( cnt_high <= M_MAX )     

     NUM = LANH0_MPI(MYID, ROOTID, N, X, VJ, CNT_LOW, CNT_HIGH, HX, ALPHA, BETA)    

     if (NUM <= 0) then
         LAN_CONV_MPI = -(CNT_LOW-1)
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
                print *, 'eig:', eig
            end if
        end if

     end if

     ALPHA(1:CNT_HIGH) = ALPHA_OLD(1:CNT_HIGH)
     BETA (1:CNT_HIGH) = BETA_OLD (1:CNT_HIGH)

     call MPI_BCAST(INFO, 2, MPI_INTEGER, ROOTID, MPI_COMM_WORLD, IERR)
      
     if ((INFO(1)==0) .and. (INFO(2)==0)) then
         call MPI_BCAST(RES, 1, MPI_DOUBLE_PRECISION, ROOTID, MPI_COMM_WORLD, IERR)
     else
         LAN_CONV_MPI = -num
         return 
     end if

     if ( RES < ETOL ) then
          LAN_CONV_MPI = cnt_high
          exit
     end if    

     LAN_CONV_MPI = - cnt_high 

     cnt_low  = cnt_high + 1
     if (cnt_low >= M_MAX) then
           lan_conv_mpi = - M_MAX
           exit
      end if
      cnt_high = cnt_high + stepCnt
      if (cnt_high > M_MAX) then
           cnt_high = M_MAX
      end if
  end do

  if (myID==rootID)  call Reorder('A', M, eig)
  
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
integer function LAN_MPI(MYID, ROOTID, N, X, M, HX, EIG)
  implicit none
  include 'mpif.h'
  integer, intent(IN)  :: MYID, ROOTID
  integer, intent(IN)  :: N, M
  double precision, intent(in) :: X(N)
  external             :: HX
  double precision, intent(out) :: EIG(M)
 

!cccccccccccccccccccc Local variables
  double precision :: BETA(M)
  integer :: LANHIJ_MPI, IERR
  integer, parameter :: LDZ = 1  
  double precision, dimension(1) :: WORK, Z    
 
  LAN_MPI = LANHIJ_MPI(MYID, ROOTID, N, X, M, HX, EIG, BETA)

  if ( LAN_MPI > 0 ) then
     if (myid == rootID) then
          call DSTEV('N', M, EIG, BETA, Z, LDZ, WORK, LAN_MPI)
     end if

     if (LAN_MPI == 0) then      
         LAN_MPI = M
     else
         LAN_MPI = -M     ! Lapack error
     end if
  end if

  call MPI_BCAST(LAN_MPI, 1, MPI_INTEGER, ROOTID, MPI_COMM_WORLD, IERR)

end 

!**************************************************************************


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    LANCZOS algorithm  to get the eigen values near energy E     c
!c       For two sets of, Apply for Symmetric matrix H             c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
integer function LANEIG_MPI(MYID, ROOTID,  N, X, M1, M2, HX, EIG1, EIG2)
  implicit none
  include 'mpif.h'
  integer, intent(IN)  :: MYID, ROOTID
  integer, intent(IN)  :: N, M1, M2
  double precision, intent(in) :: X(N)
  external             :: HX
  double precision, intent(out) :: EIG1(M1)
  double precision, intent(out) :: EIG2(M2)


!ccccccccccccccc  Local Variables cccccccccccccccccccccccccc  
  integer  :: M_MAX, INFO, IERR
  double precision, allocatable :: ALPHA(:), BETA(:),TMPB(:)
  integer, parameter :: LDZ = 1
  double precision, dimension(LDZ) :: WORK, Z
  integer :: LANHIJ_MPI

  M_MAX = max(M1, M2);   LANEIG_MPI = 0

  allocate( ALPHA(M_MAX),BETA(M_MAX), TMPB(M_MAX), STAT=INFO)
  if (INFO /= 0)    return  

  LANEIG_MPI = LANHIJ_MPI(MYID, ROOTID, N, X, M_MAX, HX, ALPHA, BETA)

  if ( LANEIG_MPI > 0 ) then
    if (MYID==ROOTID) then
      EIG1(1:M1) = ALPHA(1:M1)
      TMPB(1:M1) = BETA (1:M1)
      call DSTEV('N', M1, EIG1, TMPB, Z, LDZ, WORK, INFO)      

      if (INFO /= 0) then
           LANEIG_MPI = - M_MAX
      else
           EIG2(1:M2) = ALPHA(1:M2)
           TMPB(1:M2) = BETA(1:M2)
           call DSTEV('N',M2, EIG2, TMPB, Z, LDZ, WORK, INFO)             
           if (INFO /= 0)   LANEIG_MPI = - M_MAX           
      end if
    end if
  end if

  call MPI_BCAST(LANEIG_MPI, 1, MPI_INTEGER, ROOTID, MPI_COMM_WORLD, IERR)

  deallocate(ALPHA, BETA, TMPB, STAT=INFO)

end 
!**************************************************************************

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A Slow version but easy to implement                       c
!c         It will reuse the previous generated VJ, and Hij               c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function LAN_CONVERG_MPI(MYID, ROOTID, E0, ETOL, nType, N, X,  &
                   startCnt, stepCnt, stepEig, M, M_MAX, HX, EIG, RES)
  implicit none
  include 'mpif.h'
  integer, intent(IN)          :: MYID, ROOTID
  double precision, intent(in) :: E0,ETOL
  integer, intent(IN)          :: nType, N, M, M_MAX
  integer, intent(IN)          :: startCnt, stepCnt, stepEig
  double precision, intent(in)  :: X(N)
  double precision, intent(out) :: EIG(M)   
  external                      :: HX    
  double precision, intent(out) :: RES  


!cccccccccccccccccccccccccccccccc
  double precision :: tmp1Eig(M_MAX),tmp2Eig(M_MAX),oldEig(M)
  integer          :: cnt1, cnt2, num, IERR

!cccccccccccccccccccccccccccccccc
  integer           :: LANEIG_MPI
  double precision  :: getDiffMax, getSepMax
   
!************************************   
  LAN_CONVERG_MPI = 0
  if ((stepEig >= M_MAX) .or. (startCnt > M_MAX) .or. (M > M_MAX))  return

  Cnt1  = max(startCnt, M, stepEig+1) 

  LAN_CONVERG_MPI = - CNT1
  do  while ( CNT1 <=M_MAX )

      CNT2    = CNT1 - stepEig

      num = LANEIG_MPI( MYID, ROOTID, N, X, CNT1, CNT2,   &
                    HX, tmp1EIG, tmp2Eig)  

      if ( num <= 0) then
          LAN_CONVERG_MPI = - CNT1
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
          print *, 'eig:',eig
     end if
 
     call MPI_BCAST(RES,1,MPI_DOUBLE_PRECISION,rootID,MPI_COMM_WORLD,IERR)

     if ( RES < ETOL ) then
          LAN_CONVERG_MPI = CNT1
          exit
     end if    

     CNT1 = CNT1 + stepCnt
  end do

  if (myID == rootID)       call Reorder('A', M, eig)  
 
end 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

