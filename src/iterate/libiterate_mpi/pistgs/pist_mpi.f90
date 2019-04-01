!*****************************************************************
!*             REAL VERSION OF PIST IN MPI ENVIRONMENT           *
!*****************************************************************

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A fast version but difficult to implement                  c
!c         It will reuse the previous generated VJ, and Hij               c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function PIST_CONV_MPI(MYID, ROOTID, E0, ETOL,NTYPE, N, X,  &
          STARTCNT, STEPCNT, STEPEIG, M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  include 'mpif.h'
  integer, intent(IN)  :: MYID, ROOTID
  double precision, intent(in)  :: E0, ETOL
  integer, intent(IN        )   :: nType, N, M, M_MAX
  double precision, intent(in)  ::  X(N)  
  integer, intent(IN)           :: startCNT, stepCNT, stepEIG
  external                      :: H0X
  integer, external             :: LinSolv

  double precision, intent(out) ::EIG(M)   
  double precision, intent(out) :: RES  


!cccccccccccccccccccccccccccccccc
  double precision, dimension(N, M_MAX)     :: VJ
  double precision, dimension(M_MAX, M_MAX) :: HMAT, HMAT_OLD
  double precision, dimension(M_MAX) :: tmpEig
  double precision, dimension(M)     :: oldEig  
  
  integer   :: cnt_high, cnt_low, cnt_start
  integer   :: I, J, NUM,IERR
  integer, dimension(2) :: CONVINT 

!cccccccccccccccccccccccccccccccc  
  double precision  :: getDiffMax, getSepMax
  integer           :: PISTH0_MPI

!ccccccccccccc   Parameters to call SSTEV subroutine cccccccccccccc
  integer          :: LWORK    ! LWORK = 3*N
  double precision :: WORK(3*M_MAX)     ! 
  LWORK = 3*M_MAX  

!************************************ 
 if ((stepEig >= M_MAX) .or. (startCnt > M_MAX) .or. (M > M_MAX)) then
      PIST_CONV_MPI = 0
      return
  end if

  cnt_low   = 1
  cnt_high  = max(startCnt, M, stepEig+1) 

  CONVINT(1:2) = 0 

  do  while ( cnt_high <=M_MAX )   

     NUM = PISTH0_MPI(MYID, ROOTID, N, X,     &
               M_MAX, CNT_LOW, CNT_HIGH, H0X, LinSolv, VJ, HMAT_OLD)

     if (NUM <= 0) then
          PIST_CONV_MPI = -CNT_HIGH
          return
     end if

     if (MYID == ROOTID) then
          NUM = CNT_HIGH-STEPEIG
          HMAT(1:NUM,1:NUM) = HMAT_OLD(1:NUM,1:NUM)     
          call DSYEV('N', 'U', NUM, HMAT, M_MAX, tmpEIG, WORK,    &
                     LWORK, PIST_Conv_MPI)
         if (PIST_CONV_MPI == 0) then      
              call getWindow(E0, NUM, tmpEig, M, OLDEig) 

              NUM = cnt_high
              HMAT(1:NUM,1:NUM) = HMAT_old(1:NUM,1:NUM) 
              call DSYEV('N', 'U', NUM, HMAT, M_MAX, tmpEIG,     &
                       WORK, LWORK, PIST_Conv_MPI)
 
              if (PIST_CONV_MPI == 0) then  
                   call getWindow(E0, NUM, tmpEig, M, Eig)

                   select case (nType)
                   case (:0)
                       RES = getDiffMax(M,eig, oldEig)        

                   case (1:)
                       RES = getSepMax(M, oldEig, cnt_high, tmpEig)
                   end select
      
                   print *, 'CNT=',CNT_HIGH, ' RES=', RES

                   if ( RES <= ETOL ) then
                        CONVINT(1) = 1  
                        CONVINT(2) = cnt_high                    
                   end if  
              else                     !  Second Lapack error   
                  CONVINT(1) = 1
                  CONVINT(2) = -(cnt_HIGH - STEPEIG)
              end if
         else                 ! First Lapack error 
              CONVINT(1) = 1  
              CONVINT(2) = -cnt_high             
         end if
     end if   ! ROOTID FINISH ITS WORK

     call MPI_BCAST(CONVINT, 2, MPI_INTEGER, ROOTID, MPI_COMM_WORLD, IERR)

     PIST_CONV_MPI=CONVINT(2)
     if (CONVINT(1) == 1)    exit

     cnt_low  = cnt_high + 1
     if (cnt_low >= M_MAX) then
           pist_conv_MPI = - M_MAX
           exit
      end if
      cnt_high = cnt_high + STEPCNT
      if (cnt_high > M_MAX) then
           cnt_high = M_MAX
      end if
  end do
   
  if (MYID == ROOTID) then
      call Reorder('A', M, eig)
  end if
  
end function



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A slow version but easy to implement                       c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function PIST_CONVERG_MPI(MYID, ROOTID, E0, ETOL,NTYPE, N, X,  &
          STARTCNT, STEPCNT, STEPEIG, M, M_MAX, H0X, LinSolv, EIG, RES)
  implicit none
  include 'mpif.h'
  integer, intent(IN)           :: MYID, ROOTID
  double precision, intent(in)  :: E0, ETOL
  integer, intent(IN        )   :: nType, N, M, M_MAX
  double precision, intent(in)  ::  X(N)  
  integer, intent(IN)           :: startCNT, stepCNT, stepEIG
  external                      :: H0X
  integer, external             :: LinSolv

  double precision, intent(out) ::EIG(M)   
  double precision, intent(out) :: RES  


!cccccccccccccccccccccccccccccccc
  double precision, dimension(M_MAX) :: tmp1Eig
  double precision, dimension(M_MAX) :: tmp2Eig
  double precision, dimension(M)     :: oldEig  
  integer          ::  cnt1, cnt2, IERR

!cccccccccccccccccccccccccccccccc
  integer           :: PIST_EIG_MPI
  double precision  :: getDiffMax, getSepMax
   
!************************************   
  if ((stepEig >= M_MAX) .or. (startCnt > M_MAX) .or. (M > M_MAX)) then
      PIST_CONVERG_MPI = 0
      return
  end if

  Cnt1  = max(startCnt, M, stepEig+1) 

  do  while ( CNT1 <=M_MAX )
      CNT2    = CNT1-STEPEIG

      PIST_CONVERG_MPI = PIST_EIG_MPI( MYID, ROOTID, N, X,   &
                CNT1, CNT2, H0X, LinSolv, tmp1EIG, tmp2Eig)  

      if (PIST_CONVERG_MPI <= 0) then
          PIST_CONVERG_MPI = -CNT1
          return
      end if
       
      if (MYID == ROOTID) then
          call getWindow(E0, CNT1, tmp1Eig, M, Eig) 
          call getWindow(E0, CNT2, tmp2Eig, M, oldEig)     

          select case (nType)
          case (:0)
              RES = getDiffMax(M, eig, oldEig)        

          case (1:)
              RES = getSepMax(M, oldEig, cnt1, tmp1Eig)
          end select  

          print *, 'CNT1=', CNT1,' CNT2=', CNT2, ' RES=', RES
      end if

      call MPI_BCAST(RES, 1, MPI_DOUBLE_PRECISION, ROOTID, MPI_COMM_WORLD, IERR)
  
      if ( RES < ETOL ) then
          PIST_CONVERG_MPI = CNT1
          exit
      end if    

      CNT1 = CNT1 + STEPCNT
  end do

  if (MYID == ROOTID) then
     call Reorder('A', M, eig)
  end if
  
end function


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                 c
!c    PIST algorithm  to get the eigen values near energy E        c
!c       Apply for Symmetric matrix H                              c
!c                                                                 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
integer function PIST_MPI(MYID, ROOTID, N, X, M, H0X, LinSolv, EIG)
  implicit none
  include 'mpif.h'
  integer, intent(IN)  :: MYID, ROOTID
  integer, intent(IN)  :: N, M
  double precision, intent(in) :: X(N) 
  external             :: H0X
  integer, external    :: LinSolv

  double precision, intent(out) :: EIG(M)


!cccccccccccccccccccc Local variables
  double precision :: HMAT(M,M)   ! HMAT = <U|H|U>
  integer :: PISTHIJ_MPI, INFO, IERR
  integer          :: LWORK   ! LWORK = 3*N
  double precision :: WORK(3*M)    
 
  LWORK = 3*M     

  PIST_MPI = PISTHIJ_MPI(MYID, ROOTID, N, X, M, H0X, LinSolv, HMAT)

  if ( (PIST_MPI > 0) .and. (MYID==ROOTID) ) then
      call DSYEV('N', 'U', M, HMAT, M, EIG, WORK, LWORK, INFO)
  end if

  call MPI_BCAST(INFO, 1, MPI_DOUBLE_PRECISION, ROOTID, MPI_COMM_WORLD, IERR)

  if (INFO == 0) then      
         PIST_MPI = M
  else
         PIST_MPI = -M     ! Lapack error
  end if  

end function

!**************************************************************************


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                                 c
!c    PIST algorithm  to get the eigen values near energy E        c
!c       For two sets of 
!c       Apply for Symmetric matrix H                              c
!c                                                                 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
integer function PIST_EIG_MPI(MYID, ROOTID, N, X,     &
           M1, M2, H0X, LinSolv, EIG1, EIG2)
  implicit none
  include 'mpif.h'
  integer, intent(IN)  :: MYID, ROOTID
  integer, intent(IN)  :: N, M1, M2
  double precision, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv

  double precision, intent(out) :: EIG1(M1)
  double precision, intent(out) :: EIG2(M2)


!ccccccccccccccc  Local Variables cccccccccccccccccccccccccc  
  integer  :: M_MAX, INFO, LWORK, IERR
  double precision, allocatable :: HMAT(:,:), HOLD(:,:)  
  double precision, allocatable :: WORK(:)     
  integer :: PISTHIJ_MPI

  M_MAX = max(M1, M2)
  LWORK = 3*M_MAX  
  PIST_EIG_MPI = -1

  allocate( HMAT(M_MAX, M_MAX),HOLD(M_MAX, M_MAX),        &
            WORK(LWORK), STAT=INFO)
  if (INFO /= 0) then
     return 
  end if

  PIST_EIG_MPI = PISTHIJ_MPI(MYID, ROOTID,N,X,M_MAX,H0X,LinSolv,HMAT)

  if ( PIST_EIG_MPI > 0 .and. (MYID == ROOTID)) then
      HOLD(1:M_MAX, 1:M_MAX) = HMAT(1:M_MAX, 1:M_MAX)
      call DSYEV('N', 'U', M1, HMAT, M_MAX, EIG1, WORK, LWORK, INFO)

      if (INFO /= 0) then
           PIST_EIG_MPI = - M_MAX
      else
           call DSYEV('N', 'U', M2, HOLD, M_MAX, EIG2, WORK, LWORK, INFO)             
           if (INFO /= 0) then
                  PIST_EIG_MPI = - M_MAX
           end if
      end if
  end if

  call MPI_BCAST(PIST_EIG_MPI, 1, MPI_INTEGER, ROOTID, MPI_COMM_WORLD, IERR)

  deallocate( HOLD, HMAT, WORK, STAT=INFO)

end function

!**************************************************************************

