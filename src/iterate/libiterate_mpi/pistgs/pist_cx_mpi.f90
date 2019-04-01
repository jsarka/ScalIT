!*****************************************************************
!*     COMPLEX VERSION OF PIST IN MPI ENVIRONMENT                *
!*****************************************************************

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A fast version but difficult to implement                  c
!c         It will reuse the previous generated VJ, and Hij               c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function PIST_CONV_CX_MPI(MYID, ROOTID, E0, ETOL, NTYPE, N,  X,  &
                 startCNT,STEPcnt,stepeig,M,M_MAX,H0XCX,LinSolvCX,EIG,RES)
  implicit none
  include 'mpif.h'
  integer, intent(IN)          :: MYID, ROOTID
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, NTYPE, M, M_MAX
  double complex, intent(in) :: X(N)
  integer, intent(IN)  :: STARTCNT, STEPCNT, STEPEIG        
  external             :: H0XCX
  integer, external    :: LinSolvCX     
  double complex, intent(out) ::EIG(M)
  double precision, intent(out) :: RES              

  integer, parameter   :: nDiffTYPE    = 1  ! Real part
  integer, parameter   :: nConvType    = 6  ! Real^2+Img^2


!cccccccccccccccccccccccccccccccc
  double complex :: VJ (N, M_MAX) 
  double complex :: HMAT(M_MAX, M_MAX), HMAT_OLD(M_MAX,M_MAX)
  double complex :: tmpEig(M_MAX),oldEig(M)  
  
  integer   :: cnt_high, cnt_low 
  integer   :: I, J, NUM,IERR, CONVINT(2)

!cccccccccccccccccccccccccccccccc  
  double precision  :: getDiffMax_CX, getSepMAX_CX
  integer           :: PISTH0_CX_MPI

!ccccccccccccc   Parameters to call ZGEEV subroutine cccccccccccccc
  integer, parameter ::  LDVL  = 1,LDVR = 1, SCALE = 3
  double complex     ::  VL(LDVL), VR(LDVR)
  integer            ::  LWORK 
  double complex  :: Work(SCALE*M),RWork(scale*N)
    
  LWORK = SCALE*M_MAX  


!************************************ 
  if ((stepEig >= M_MAX) .or. (startCnt > M_MAX) .or. (M > M_MAX)) then
      PIST_CONV_CX_MPI = 0
      return
  end if

  cnt_low   = 1
  cnt_high  = max(startCnt, M, stepEig+1)

  CONVINT(1:2) = 0 

  do  while ( cnt_high <=M_MAX )
   
     NUM = PISTH0_CX_MPI(MYID, ROOTID, N, X,  M_MAX,  &
                CNT_LOW, CNT_HIGH, H0XCX, LinSolvCX, VJ, HMAT_OLD)

     if (NUM <= 0) then
          PIST_CONV_CX_MPI = -CNT_HIGH
          return
     end if

     if (MYID == ROOTID) then
          num = cnt_high-stepEig
          HMAT(1:num,1:num) = HMAT_OLD(1:num,1:num) 

          call ZGEEV('N', 'N', NUM, HMAT, M_MAX, tmpEIG, VL, LDVL, VR, LDVR, &
                     WORK, LWORK, RWORK, PIST_Conv_CX_MPI)

         if (PIST_CONV_CX_MPI == 0) then      

              call getWindow_CX(E0,nDiffType, num, tmpEig, M, oldEig) 

              num = cnt_high
              HMAT(1:num,1:num) = HMAT_OLD(1:num,1:num)      

              call ZGEEV('N', 'N', NUM, HMAT, M_MAX, tmpEIG, VL, LDVL, VR, LDVR, &
                       WORK, LWORK, RWORK, PIST_Conv_CX_MPI)
 
              if (PIST_CONV_CX_MPI == 0) then  

                   call getWindow_CX(E0, nDiffType, num, tmpEig, M, Eig)

                   select case (nType)
                   case (:0)
                        RES = getDiffMax_CX(nConvTYPE, M,eig, oldEig)        

                   case (1:)
                        RES = getSepMax_CX(nConvTYPE, M, oldEig, cnt_high, tmpEig)
                   end select
      
                   print *, 'CNT1=',CNT_HIGH, ' Cnt2=',(cnt_High-stepEig),' RES=', RES
!                   print *, 'eig:',oldEig

                   if ( RES <= ETOL ) then
                        CONVINT(1) = 1  
                        CONVINT(2) = cnt_high                    
                   end if  
              else                     !  Second Lapack error   
                  CONVINT(1) = 1
                  CONVINT(2) = -(cnt_high-stepEig)          
              end if
         else                 ! First Lapack error 
              CONVINT(1) = 1  
              CONVINT(2) = -cnt_high             
         end if
     end if   ! ROOTID FINISH ITS WORK

     call MPI_BCAST(CONVINT, 2, MPI_INTEGER, ROOTID, MPI_COMM_WORLD, IERR)

     PIST_CONV_CX_MPI=CONVINT(2)

     if (CONVINT(1) == 1)   exit

     cnt_low  = cnt_high + 1
     if (cnt_low >= M_MAX) then
           pist_conv_CX_MPI = - M_MAX
           exit
      end if
      cnt_high = cnt_high + STEPCNT

      if (cnt_high > M_MAX)       cnt_high = M_MAX   

  end do
 
  if (MYID == ROOTID) then
      call Reorder_CX('A', NDIFFTYPE,  M, eig)
  end if

end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A slow version but easy to implement                       c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function PIST_CONVERG_CX_MPI(MYID, ROOTID, E0, ETOL, NTYPE, N,  X, &
                 startCNT,STEPcnt,stepeig,M,M_MAX,H0XCX,LinSolvCX,EIG,RES)
  implicit none
  include 'mpif.h'
  integer, intent(IN)          :: MYID, ROOTID
  double precision, intent(in) :: E0, ETOL
  integer, intent(IN)  :: N, NTYPE, M, M_MAX
  double complex, intent(in) :: X(N)
  integer, intent(IN)  :: STARTCNT, STEPCNT, STEPEIG        
  external             :: H0XCX
  integer, external    :: LinSolvCX

  double complex, intent(out) ::EIG(M)
  double precision, intent(out) :: RES              

  integer, parameter :: nDiffTYPE    = 1  ! Real part
  integer, parameter :: nConvType    = 6  ! Real^2+Img^2


!cccccccccccccccccccccccccccccccc
  double complex :: tmp1Eig(M_MAX), tmp2Eig(M_MAX), oldEig(M)
  integer        ::  cnt1, cnt2, IERR

!cccccccccccccccccccccccccccccccc
  integer           :: PIST_EIG_CX_MPI
  double precision  :: getDiffMax_CX, GETSEPMAX_CX
   
!************************************   
  if ((stepEig >= M_MAX) .or. (startCnt > M_MAX) .or. (M > M_MAX)) then
      PIST_CONVERG_CX_MPI = 0
      return
  end if

  Cnt1  = max(startCnt, M, stepEig+1) 

  do  while ( CNT1 <=M_MAX )
      CNT2    = CNT1-STEPEIG

      PIST_CONVERG_CX_MPI = PIST_EIG_CX_MPI(MYID, ROOTID, N, X,    &
             CNT1, CNT2, H0XCX, LinSolvCX,tmp1EIG,tmp2Eig)  

      if (PIST_CONVERG_CX_MPI <= 0) then
          PIST_CONVERG_CX_MPI = - CNT1
          return
      end if
       
      if (MYID == ROOTID) then
          call getWindow_CX(E0, nType, CNT1, tmp1Eig, M, Eig) 
          call getWindow_CX(E0, nType, CNT2, tmp2Eig, M, oldEig)     
      
          select case (nType)
          case (:0)
               RES = getDiffMax_CX(nConvType,M, eig, oldEig)        

          case (1:)
               RES = getSepMax_CX(nConvTYPE, M, oldEig, cnt1, tmp1Eig)
          end select  

          print *, 'CNT1=',CNT1,' Cnt2=',cnt2,' RES=', RES
      end if

      call MPI_BCAST(RES, 1, MPI_DOUBLE_PRECISION, ROOTID, MPI_COMM_WORLD, IERR)
  
      if ( RES < ETOL ) then
          PIST_CONVERG_CX_MPI = CNT1
          exit
      end if    

      CNT1 = CNT1 + STEPCNT
  end do
 
  if (MYID == ROOTID)  call Reorder_CX('A', NDIFFTYPE, M, eig) 
  
end 


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    PIST algorithm  to get the eigen values near energy E        c
!c       Apply for Symmetric matrix H                              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
integer function PIST_CX_MPI(MYID, ROOTID, N, X, M, H0X, LinSolv, EIG)
  implicit none
  include 'mpif.h'
  integer, intent(IN)  :: MYID, ROOTID
  integer, intent(IN)  :: N, M                       
  double complex, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv

  double complex, intent(out) :: EIG(M)


!cccccccccccccccccccc Local variables
  double complex, dimension(M, M) :: HMAT   ! HMAT = <U|H|U>
  integer :: PISTHIJ_CX_MPI, INFO, IERR
  integer, parameter   ::  LDVL = 1, LDVR = 1,SCALE = 3
  double complex :: VL(LDVL),VR(LDVR), WORK(SCALE*M)
  integer        ::  LWORK 
  double precision  ::  RWORK(2*N)
 
  LWORK = SCALE*M     

  PIST_CX_MPI = PISTHIJ_CX_MPI(MYID, ROOTID, N, X, M, H0X, LinSolv, HMAT)

  if ( (PIST_CX_MPI > 0) .and. (MYID==ROOTID) ) then
      call ZGEEV('N', 'N', M, HMAT, M, EIG, VL, LDVL, VR, LDVR, WORK,      &
                   LWORK, RWORK, INFO)
  end if

  call MPI_BCAST(INFO, 1, MPI_DOUBLE_PRECISION, ROOTID, MPI_COMM_WORLD, IERR)

  if (INFO == 0) then      
         PIST_CX_MPI = M
  else
         PIST_CX_MPI = -M     ! Lapack error
  end if  

end 

!**************************************************************************


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    PIST algorithm  to get the eigen values near energy E        c
!c       For two sets,   Apply for Symmetric matrix H              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
integer function PIST_EIG_CX_MPI(MYID,ROOTID,N,X,M1,M2,H0X,LinSolv,EIG1,EIG2)
  implicit none
  include 'mpif.h'
  integer, intent(IN)  :: MYID, ROOTID
  integer, intent(IN)  :: N, M1, M2
  double complex, intent(in) :: X(N)
  external             :: H0X
  integer, external    :: LinSolv

  double complex, intent(out) :: EIG1(M1), EIG2(M2)


!ccccccccccccccc  Local Variables cccccccccccccccccccccccccc  
  integer  :: M_MAX, INFO, LWORK, IERR
  double complex, allocatable :: HMAT(:,:), HOLD(:,:), WORK(:)
  double precision, allocatable  :: RWORK(:)
  integer :: PISTHIJ_CX_MPI

  integer, parameter  ::  LDVL = 1, LDVR = 1,SCALE = 3 
  double complex :: VL(LDVL), VR(LDVR)

  M_MAX = max(M1, M2);    LWORK = 3*M_MAX  
  PIST_EIG_CX_MPI = -1

  allocate(HMAT(M_MAX, M_MAX),HOLD(M_MAX, M_MAX),        &
            WORK(LWORK), RWORK(2*M_MAX), STAT=INFO)
  if (INFO /= 0)     return   

  PIST_EIG_CX_MPI = PISTHIJ_CX_MPI(MYID,ROOTID,N,X,M_MAX,H0X,LinSolv,HMAT)

  if ( PIST_EIG_CX_MPI > 0 .and. (MYID == ROOTID)) then
      HOLD(1:M_MAX, 1:M_MAX) = HMAT(1:M_MAX, 1:M_MAX)
      call ZGEEV('N', 'N', M1, HMAT, M_MAX, EIG1, VL, LDVL, VR, LDVR, WORK,      &
                   LWORK, RWORK, INFO)

      if (INFO /= 0) then
           PIST_EIG_CX_MPI = - M1
      else
          call ZGEEV('N', 'N', M2, HOLD, M_MAX, EIG2, VL, LDVL, VR, LDVR, WORK, &
                   LWORK, RWORK, INFO)
           if (INFO /= 0)     PIST_EIG_CX_MPI = - M2           
      end if
  end if

  call MPI_BCAST(PIST_EIG_CX_MPI, 1, MPI_INTEGER, ROOTID, MPI_COMM_WORLD, IERR)

  deallocate( HOLD, HMAT, WORK, RWORK, STAT=INFO)

end 

!**************************************************************************

