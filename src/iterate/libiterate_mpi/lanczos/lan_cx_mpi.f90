!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A fast version but difficult to implement                  c
!c         It will reuse the previous generated VJ, and Hij               c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function lanConv_CX_MPI(nodNum, myID, rootID, E0, Etol, nType, &
                                N, X, startCnt, stepCnt, stepEig, M,   &
                                M_MAX, HX, EIG, RES)
  implicit none
  include 'mpif.h'
  integer, intent(IN)          :: nodNum, myID, rootID
  double precision, intent(in) :: E0,Etol
  integer, intent(IN)          :: nType, N, M, M_MAX
  integer, intent(IN)          :: startCnt, stepCnt, stepEig
  double precision,intent(in)  :: X(N)  
  double precision,intent(out) :: EIG(M)   
  external                     :: HX    
  double precision, intent(out) :: RES             


!cccccccccccccccccccccccccccccccc
  double precision :: VJ(N, M_MAX+1)
  double precision, dimension(M_MAX) :: ALPHA, ALPHA_OLD
  double precision, dimension(M_MAX) :: BETA, BETA_OLD
  double precision :: oldEig(M)  
  
  integer   :: cnt_high, cnt_low, cnt_start
  integer   :: I, J, NUM, INFO(2)

!cccccccccccccccccccccccccccccccc  
  double precision  :: getDiffMax, getSepMax
  integer           :: LANH0_MPI
  double precision  :: tmp

!ccccccccccccc   Parameters to call SSTEV subroutine cccccccccccccc
  integer, parameter   :: LDZ = 1   
  double precision, dimension(LDZ) :: WORK, Z

!************************************ 
  lanConv_CX_MPI = 0
  if ((stepEig >= M_MAX) .or. (startCnt > M_MAX) .or. (M > M_MAX)) &
      return

  ALPHA(1:M_MAX) = 0.0D0
  BETA(1:M_MAX)  = 0.0D0

  cnt_low   = 1
  Cnt_high  = max(startCnt, M, stepEig+1)   

  lanConv_CX_MPI = -Cnt_high

  do  while ( cnt_high <= M_MAX )     

     NUM = LANH0_MPI(nodNum, myID, rootID, N, X, VJ, CNT_LOW, CNT_HIGH, HX, ALPHA, BETA)
     if (NUM <= 0) then
         lanConv_CX_MPI = -(CNT_LOW-1)
         return
     end if

     if (myID == rootID) then
        ALPHA_OLD(1:CNT_HIGH) = ALPHA(1:CNT_HIGH)
        BETA_OLD (1:CNT_HIGH) = BETA (1:CNT_HIGH)
    
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
            end if
        end if
     end if

     call MPI_BCAST(INFO, 2, MPI_INTEGER, rootID, MPI_COMM_WORLD)
      
     if ((INFO(1)==0) .and. (INFO(2)==0)) then
         call MPI_BCAST(RES, 1, MPI_DOUBLE_PRECISION, rootID, MPI_COMM_WORLD)
     else
         lanConv_CX_MPI = -num
         return 
     end if

     if ( RES < Etol ) then
          lanConv_CX_MPI = cnt_high
          exit
     end if    

     ALPHA(1:CNT_HIGH) = ALPHA_OLD(1:CNT_HIGH)
     BETA (1:CNT_HIGH) = BETA_OLD (1:CNT_HIGH)

     lanConv_CX_MPI = - cnt_high 

     cnt_low  = cnt_high + 1
     if (cnt_low >= M_MAX) then
           lanConv_CX_MPI = - M_MAX
           exit
      end if
      cnt_high = cnt_high + stepCnt
      if (cnt_high > M_MAX) cnt_high = M_MAX
  end do

  call Reorder('A', M, eig)
  
end function lanConv_CX_MPI



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                 c
!c    LANCZOS algorithm  to get the eigen values near energy E        c
!c       Apply for Symmetric matrix H                              c
!c                                                                 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Get the eigen values via Full Reorthogonization Lanczos:
!
! integer function Lan(N, X, M, HX, EIG)
!
!
integer function Lan_CX_MPI(nodNum, myID, rootID, N, X, M, HX, EIG)
  implicit none
  INCLUDE 'mpif.h'
  integer, intent(IN)  :: nodNum, myID, rootID
  integer, intent(IN)  :: N, M
  double precision, intent(in) :: X(N)
  external             :: HX

  double precision, intent(out) :: EIG(M)
 

!cccccccccccccccccccc Local variables
  double precision :: BETA(M)  
  integer :: LANHIJ_MPI
  integer, parameter :: LDZ = 1  
  double precision   :: WORK(LDZ), Z(LDZ)    
 
  Lan_CX_MPI = LANHIJ_MPI(nodNum, myID, rootID, N, X, M, HX, EIG, BETA)

  if ( Lan_CX_MPI > 0 ) then
     if (myid == rootID) then
          call DSTEV('N', M, EIG, BETA, Z, LDZ, WORK, Lan_CX_MPI)
     end if

     if (Lan_CX_MPI == 0) then      
         Lan_CX_MPI = M
     else
         Lan_CX_MPI = -M     ! Lapack error
     end if
  end if

end function Lan_CX_MPI

!**************************************************************************


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                 c
!c    LANCZOS algorithm  to get the eigen values near energy E        c
!c       For two sets of 
!c       Apply for Symmetric matrix H                              c
!c                                                                 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
integer function lanEig_CX_MPI(nodNum, myID, rootID,  N, X, M1, M2, HX, EIG1, EIG2)
  implicit none
  integer, intent(IN)  :: nodNum, myID, rootID
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
  integer :: LANHIJ_MPI

  M_MAX = max(M1, M2)
  lanEig_CX_MPI = 0

  allocate( ALPHA(M_MAX),BETA(M_MAX), TMPB(M_MAX), STAT=INFO)
  if (INFO /= 0)   return   

  lanEig_CX_MPI = LANHIJ_MPI(nodNum, myID, rootID, N, X, M_MAX, HX, ALPHA, BETA)

  if ( lanEig_CX_MPI > 0 ) then
    if (myID==rootID) then
      EIG1(1:M1) = ALPHA(1:M1)
      TMPB(1:M1) = BETA (1:M1)
      call DSTEV('N', M1, EIG1, TMPB, Z, LDZ, WORK, INFO)      

      if (INFO /= 0) then
           lanEig_CX_MPI = - M_MAX
      else
           EIG2(1:M2) = ALPHA(1:M2)
           TMPB(1:M2) = BETA(1:M2)
           call DSTEV('N',M2, EIG2, TMPB, Z, LDZ, WORK, INFO)             
           if (INFO /= 0) lanEig_CX_MPI = - M_MAX         
      end if
    end if
  end if

  deallocate(ALPHA, BETA, TMPB, STAT=INFO)

end function lanEig_CX_MPI

!**************************************************************************


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A Slow version but easy to implement                       c
!c         It will reuse the previous generated VJ, and Hij               c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function lan_Converg_CX_MPI(nodNum, myID, rootID, E0, Etol, &
                                    nType, N, X, startCnt, stepCnt, &
                                    stepEig, M, M_MAX, HX, EIG, RES)
  implicit none
  include 'mpif.h'
  integer, intent(IN)          :: nodNum, myID, rootID
  double precision, intent(in) :: E0,Etol
  integer, intent(IN)          :: nType, N, M, M_MAX
  integer, intent(IN)          :: startCnt, stepCnt, stepEig
  double precision, intent(in)  :: X(N)
  double precision, intent(out) :: EIG(M)   
  external                      :: HX    
  double precision, intent(out) :: RES  


!cccccccccccccccccccccccccccccccc
  double precision :: tmp1Eig(M_MAX), tmp2Eig(M_MAX),oldEig(M)  
  integer          ::  cnt1, cnt2, num, IERR

!cccccccccccccccccccccccccccccccc
  integer           :: lanEig_CX_MPI
  double precision  :: getDiffMax, getSepMax
   
!************************************   
  lan_Converg_CX_MPI = 0
  if ((stepEig >= M_MAX) .or. (startCnt > M_MAX) .or. (M > M_MAX)) &
      return  

  Cnt1  = max(startCnt, M, stepEig+1) 

  lan_Converg_CX_MPI = - CNT1
  do  while ( CNT1 <=M_MAX )

      CNT2    = CNT1 - stepEig

      num = lanEig_CX_MPI(nodNum, myID, rootID, N, X, CNT1, CNT2,   &
                    HX, tmp1EIG, tmp2Eig)  

      if ( num <= 0) then
          lan_Converg_CX_MPI = - CNT1
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

     if ( RES < Etol ) then
          lan_Converg_CX_MPI = CNT1
          exit
     end if    

     CNT1 = CNT1 + stepCnt
  end do

  if (myID == rootID) then
      call Reorder('A', M, eig)  
  end if

end function lan_Converg_CX_MPI

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Perform Lanczos algorithm and do convergent testing
!     with full Reorthogonization.
!
! LAN_CONV:   fast, but need more memory
! LAN_CONVERG slow, but may need less memory
! 
! INTEGER FUNCTION LAN_CONV   (nodNum, myID, rootID, E0, Etol, nType, N, X, & 
!                  startCnt, stepCnt, stepEig, M, M_MAX, HX, EIG, RES)
! INTEGER FUNCTION LAN_CONVERG(nodNum, myID, rootID, E0, Etol, nType, N, X, &
!                  startCnt, stepCnt, stepEig, M, M_MAX, HX, EIG, RES)
!
! Input parameters:
!    nodNum: number of nodes
!    myID:   my id in the MPI WORLD
!    rootID: The ID of the root machine
!    E0:  Central of energy to interested in
!    Etol:Convergenc Criteria
!    nType: convergence Type, current only 0, and 1
!          <= 0 : 1 to 1 testing
!          >= 1 : 1 to many testing
!    N:   Dimension of vector
!    X:   Initial Vector
!    startCnt: start number of Lanczos Iteration
!    stepCnt:  step of Lanczos iteration for next iteration
!    stepEig:  The comparason between eigen value for convergence
!    M:        Number of interested eigen values
!    M_MAX:    Max. number of Lanczos iteration
!    HX:  External function to perform H*X, H is the real symmetric matrix
!
! Output parameters:
! EIG: The output eigen values around E0
! RES: The maximum error for convergence
! Lanczos_Conv: 
!             > 0, successful, the number of Lanczos iteration
!             < 0, fails, either Lanczos fails (<M_MAX) or convergence doesn't meet
!             = 0, Lapack Error
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
