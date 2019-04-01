!*****************************************************************
!*             REAL VERSION OF PIST IN MPI ENVIRONMENT           *
!*****************************************************************

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A fast version but difficult to implement                  c
!c         It will reuse the previous generated VJ, and Hij               c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function PISTF_CONV_MPI(MYID, ROOTID, E0, ETOL,NTYPE, N, X,  &
          STARTCNT, STEPCNT, STEPEIG, M, M_MAX, H0X, LinSolv, EIG,   &
          RES, fname,pos,gSize)
  implicit none
  include 'mpif.h'
  integer, intent(IN)  :: MYID, ROOTID
  double precision, intent(in)  :: E0, ETOL
  integer, intent(IN        )   :: nType, N, M, M_MAX
  double precision, intent(in)  ::  X(N)  
  integer, intent(IN)           :: startCNT, stepCNT, stepEIG
  external                      :: H0X
  integer, external             :: LinSolv

  double precision, intent(out) :: EIG(M)   
  double precision, intent(out) :: RES 
  character(len=*), intent(IN)  :: fname 
  integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos, gSize

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

  integer :: dbSize, intSize, fh
  integer(kind=MPI_OFFSET_KIND) :: offset


!ccccccccccccc   Parameters to call SSTEV subroutine cccccccccccccc
  integer          :: LWORK    ! LWORK = 3*N
  double precision :: WORK(3*M_MAX)     ! 
  LWORK = 3*M_MAX  

!************************************ 
 if ((stepEig >= M_MAX) .or. (startCnt > M_MAX) .or. (M > M_MAX)) then
      PISTF_CONV_MPI = 0
      return
  end if

  cnt_low   = 1
  cnt_high  = max(startCnt, M, stepEig+1) 

  CONVINT(1:2) = 0 

  do  while ( cnt_high <=M_MAX )   

     NUM = PISTH0_MPI(MYID, ROOTID, N, X,     &
               M_MAX, CNT_LOW, CNT_HIGH, H0X, LinSolv, VJ, HMAT_OLD)

     if (NUM <= 0) then
          PISTF_CONV_MPI = -CNT_HIGH
          return
     end if

     if (MYID == ROOTID) then
          NUM = CNT_HIGH-STEPEIG
          HMAT(1:NUM,1:NUM) = HMAT_OLD(1:NUM,1:NUM)     
          call DSYEV('N', 'U', NUM, HMAT, M_MAX, tmpEIG, WORK,    &
                     LWORK, PISTF_Conv_MPI)
         if (PISTF_CONV_MPI == 0) then      
              call getWindow(E0, NUM, tmpEig, M, OLDEig) 

              NUM = cnt_high
              HMAT(1:NUM,1:NUM) = HMAT_old(1:NUM,1:NUM) 
              call DSYEV('N', 'U', NUM, HMAT, M_MAX, tmpEIG,     &
                       WORK, LWORK, PISTF_Conv_MPI)
 
              if (PISTF_CONV_MPI == 0) then  
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

     PISTF_CONV_MPI=CONVINT(2)
     if (CONVINT(1) == 1)    exit

     cnt_low  = cnt_high + 1
     if (cnt_low >= M_MAX) then
           pistf_conv_MPI = - M_MAX
           exit
      end if
      cnt_high = cnt_high + STEPCNT
      if (cnt_high > M_MAX) then
           cnt_high = M_MAX
      end if
  end do

  if (PISTF_CONV_MPI>0) then  ! store data
    call MPI_TYPE_Size(MPI_DOUBLE_PRECISION, dbSize, ierr)
    call MPI_TYPE_Size(MPI_INTEGER, intSize, ierr)

    call MPI_FILE_OPEN(MPI_COMM_WORLD,fname,MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                       MPI_INFO_NULL,fh,ierr)
    if (myid==rootID) then
       offset=0;
       call MPI_File_Write_At(fh, offset, PISTF_CONV_MPI, 1, &
                   MPI_INTEGER, MPI_STATUS_IGNORE, ierr)

       offset=offset+intSize
       call MPI_File_Write_At(fh, offset, N, 1, &
                   MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
       
       offset=offset+intSize
       call MPI_File_Write_At(fh, offset, tmpEig, PISTF_CONV_MPI,&
                   MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)

       offset=offset+PISTF_CONV_MPI*dbSize
       do i = 1, PISTF_CONV_MPI
          call MPI_File_Write_At(fh, offset, HMat_old(1,i), PISTF_CONV_MPI, &
                   MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
          offset=offset+PISTF_CONV_MPI*dbSize
       end do
    endif
    offset=2*intSize+(PISTF_CONV_MPI+PISTF_CONV_MPI**2)*dbSize
    offset=offset+(pos-1)*dbSize

    do i = 1, PISTF_CONV_MPI
          call MPI_File_Write_At(fh, offset, VJ(1,i), N, &
                   MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
          offset=offset+gSize*dbSize
    end do

    call MPI_FILE_CLOSE(fh, ierr)

  end if

  if (myID==rootID) then
      if (PISTF_CONV_MPI>0) call Reorder('A', M, eig)
  end if
     
!  if (MYID == ROOTID) then
!     if (PISTF_CONV_MPI>0) then  ! store data
!        open(99, FILE=fname, status='Replace', FORM='UNFORMATTED')
!        write(99) PISTF_CONV_MPI, N
!        write(99) tmpEig(1:PISTF_CONV_MPI),     &
!                  HMat_old(1:PISTF_CONV_MPI,1:PISTF_CONV_MPI)
!        write(99) VJ(1:N, 1:PISTF_CONV_MPI)
!        close(99)
!      end if
!
!      call Reorder('A', M, eig)
!  end if
! 
end function

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!**************************************************************************

