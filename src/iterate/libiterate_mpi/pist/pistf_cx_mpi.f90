!*****************************************************************
!*     COMPLEX VERSION OF PIST IN MPI ENVIRONMENT                *
!*****************************************************************

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             A fast version but difficult to implement                  c
!c         It will reuse the previous generated VJ, and Hij               c
!c         Apply PIST until M eigen values around E0 convergence          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function PISTF_CONV_CX_MPI(MYID, ROOTID, E0, ETOL, NTYPE, N, X,  &
                 startCNT, STEPcnt, stepeig, M, M_MAX, H0XCX, LinSolvCX,    &
                 EIG, RES, fname, pos, gsize)
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
  character(LEN=*), intent(IN)  :: fname

  integer(kind=MPI_OFFSET_KIND), intent(IN) :: pos, gSize

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

  integer :: cxSize, intSize, fh
  integer(kind=MPI_OFFSET_KIND) :: offset

!ccccccccccccc   Parameters to call ZGEEV subroutine cccccccccccccc
  integer, parameter ::  LDVL  = 1,LDVR = 1, SCALE = 3
  double complex     ::  VL(LDVL), VR(LDVR)
  integer            ::  LWORK 
  double complex  :: Work(SCALE*M),RWork(scale*N)
    
  LWORK = SCALE*M_MAX  


!************************************ 
  if ((stepEig >= M_MAX) .or. (startCnt > M_MAX) .or. (M > M_MAX)) then
      PISTF_CONV_CX_MPI = 0
      return
  end if

  cnt_low   = 1
  cnt_high  = max(startCnt, M, stepEig+1)

  CONVINT(1:2) = 0 

  do  while ( cnt_high <=M_MAX )
   
     NUM = PISTH0_CX_MPI(MYID, ROOTID, N, X,  M_MAX,  &
                CNT_LOW, CNT_HIGH, H0XCX, LinSolvCX, VJ, HMAT_OLD)

     if (NUM <= 0) then
          PISTF_CONV_CX_MPI = -CNT_HIGH
          return
     end if

     if (MYID == ROOTID) then
          num = cnt_high-stepEig
          HMAT(1:num,1:num) = HMAT_OLD(1:num,1:num) 

          call ZGEEV('N', 'N', NUM, HMAT, M_MAX, tmpEIG, VL, LDVL, VR, LDVR, &
                     WORK, LWORK, RWORK, PISTF_Conv_CX_MPI)

         if (PISTF_CONV_CX_MPI == 0) then      

              call getWindow_CX(E0,nDiffType, num, tmpEig, M, oldEig) 

              num = cnt_high
              HMAT(1:num,1:num) = HMAT_OLD(1:num,1:num)      

              call ZGEEV('N', 'N', NUM, HMAT, M_MAX, tmpEIG, VL, LDVL, VR, LDVR, &
                       WORK, LWORK, RWORK, PISTF_Conv_CX_MPI)
 
              if (PISTF_CONV_CX_MPI == 0) then  

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

     PISTF_CONV_CX_MPI=CONVINT(2)

     if (CONVINT(1) == 1)   exit

     cnt_low  = cnt_high + 1
     if (cnt_low >= M_MAX) then
           pistf_conv_CX_MPI = - M_MAX
           exit
      end if
      cnt_high = cnt_high + STEPCNT

      if (cnt_high > M_MAX)       cnt_high = M_MAX   

  end do


  if (PISTF_CONV_CX_MPI>0) then  ! store data
    call MPI_TYPE_Size(MPI_DOUBLE_COMPLEX, cxSize, ierr)
    call MPI_TYPE_Size(MPI_INTEGER, intSize, ierr)

    call MPI_FILE_OPEN(MPI_COMM_WORLD,fname,MPI_MODE_WRONLY+MPI_MODE_CREATE, &
                       MPI_INFO_NULL,fh,ierr)
    if (myid==rootID) then
       offset=0;
       call MPI_File_Write_At(fh, offset, PISTF_CONV_CX_MPI, 1, &
                   MPI_INTEGER, MPI_STATUS_IGNORE, ierr)

       offset=offset+intSize
       call MPI_File_Write_At(fh, offset, N, 1, &
                   MPI_INTEGER, MPI_STATUS_IGNORE, ierr)

       offset=offset+intSize
       call MPI_File_Write_At(fh, offset, tmpEig, PISTF_CONV_CX_MPI,&
                   MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)

       offset=offset+PISTF_CONV_CX_MPI*cxSize
       do i = 1, PISTF_CONV_CX_MPI
          call MPI_File_Write_At(fh, offset, HMat_old(1,i), PISTF_CONV_CX_MPI, &
                   MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr)
          offset=offset+PISTF_CONV_CX_MPI*cxSize
       end do
    endif
    offset=2*intSize+(PISTF_CONV_CX_MPI+PISTF_CONV_CX_MPI**2)*cxSize
    offset=offset+(pos-1)*cxSize

    do i = 1, PISTF_CONV_CX_MPI
          call MPI_File_Write_At(fh, offset, VJ(1,i), N, &
                   MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr)
          offset=offset+gSize*cxSize
    end do

    call MPI_FILE_CLOSE(fh, ierr)

  end if

  if (myID==rootID) then
      if (PISTF_CONV_CX_MPI>0) call Reorder_CX('A', NDIFFTYPE,  M, eig)
  end if


 
!  if (MYID == ROOTID) then
!     if (PISTF_CONV_CX_MPI>0) then  ! store data
!        open(99, FILE=fname, status='Replace', FORM='UNFORMATTED')
!        write(99) PISTF_CONV_CX_MPI, N
!        write(99) tmpEig(1:PISTF_CONV_CX_MPI),     &
!                  HMat_old(1:PISTF_CONV_CX_MPI,1:PISTF_CONV_CX_MPI)
!        write(99) VJ(1:N, 1:PISTF_CONV_CX_MPI)
!        close(99)
!      end if
!
!      call Reorder_CX('A', NDIFFTYPE,  M, eig)
!  end if


end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
