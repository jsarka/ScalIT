!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c          Grant-Schmidt Orthogonization in MPI              c
!c       Just for one step,this is faster than MGS_MPI        c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

logical function GS_ORTH_CX_MPI(N, WJ, M, VJ, ierr)
   implicit none
   include 'mpif.h'
   integer, intent(IN) :: N, M
   double complex, intent(IN)  :: WJ(N)
   double complex, intent(INOUT) :: VJ(N,M+1)
   integer, intent(out) :: IERR


!ccccccccccccccccc   
   double complex   :: dotwv(M)
   double precision :: normw, NORM_CX_MPI
   integer :: i, M0

   M0 = M + 1
   VJ(1:N, M0) = WJ(1:N)   
   call DOTPROD_MVCX_MPI(MPI_COMM_WORLD, N, M, VJ, 1, WJ, dotwv, ierr )

   do I = 1, M
       VJ(1:N, M0) = VJ(1:N, M0) - dotwv(I) * VJ(1:N, I)
   end do

   normw = NORM_CX_MPI(MPI_COMM_WORLD, N, VJ(1:N, M0), IERR)

   if (normw == 0.0D0) then
      GS_ORTH_CX_MPI = .false.
   else
      VJ(1:N, M0) = VJ(1:N, M0) / normw
      GS_ORTH_CX_MPI = .true.
   end if

end function GS_ORTH_CX_MPI

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Modified Grant-Schimit Orthogonization           c
!c                  Just for one step                         c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

logical function MGS_ORTH_CX_MPI(N, WJ, M, VJ, IERR)
   implicit none
   include 'mpif.h'
   integer, intent(IN) :: N, M
   double complex, intent(IN)  :: WJ(N)
   double complex, intent(INOUT) :: VJ(N,M+1)
   integer, intent(out) :: IERR


!ccccccccccccccccc   
   double complex   :: dotwv, DOTPROD_SVCX_MPI
   double precision :: normw, NORM_CX_MPI
   integer :: i, M0   

   M0 = M + 1
   VJ(1:N, M0) = WJ(1:N)   
  
   do I = 1, M
      dotwv = DOTPROD_SVCX_MPI(MPI_COMM_WORLD, N, VJ(1:N, I),VJ(1:N, M0), IERR)                
      VJ(1:N, M0) = VJ(1:N, M0) - dotwv * VJ(1:N, I)
   end do
  
   normw = NORM_CX_MPI(MPI_COMM_WORLD, N, VJ(1:N, M0), IERR)

   if (normw == 0.0D0) then      
      MGS_ORTH_CX_MPI = .false.
   else
      VJ(1:N, M0) = VJ(1:N, M0) / normw
      MGS_ORTH_CX_MPI = .true.
   end if

end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

logical function LAN_GS_ORTH_CX_MPI(N, WJ, M, VJ, beta, ierr)
   implicit none
   include 'mpif.h'
   integer, intent(IN) :: N, M
   double complex, intent(IN)  :: WJ(N)
   double complex, intent(INOUT) :: VJ(N,M+1)
   integer, intent(out) :: IERR
   double precision, intent(out) :: beta


!ccccccccccccccccc   
   double complex   :: dotwv(M)
   double precision :: normw, NORM_CX_MPI
   integer :: i, M0

   M0 = M + 1
   VJ(1:N, M0) = WJ(1:N)   
   call DOTPROD_MVCX_MPI(MPI_COMM_WORLD, N, M, VJ, 1, WJ, dotwv, ierr )

   do I = 1, M
       VJ(1:N, M0) = VJ(1:N, M0) - dotwv(I) * VJ(1:N, I)
   end do

   normw = NORM_CX_MPI(MPI_COMM_WORLD, N, VJ(1:N, M0), IERR)
   beta = normw

   if (normw == 0.0D0) then
      LAN_GS_ORTH_CX_MPI = .false.
   else
      VJ(1:N, M0) = VJ(1:N, M0) / normw
      LAN_GS_ORTH_CX_MPI = .true.
   end if

end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Modified Grant-Schimit Orthogonization           c
!c                  Just for one step                         c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

logical function LAN_MGS_ORTH_CX_MPI(N, WJ, M, VJ, beta, IERR)
   implicit none
   include 'mpif.h'
   integer, intent(IN) :: N, M
   double complex, intent(IN)  :: WJ(N)
   double complex, intent(INOUT) :: VJ(N,M+1)
   double precision, intent(OUT) :: beta   
   integer, intent(out) :: IERR


!ccccccccccccccccc   
   double complex   :: dotwv, DOTPROD_SVCX_MPI
   double precision :: normw, NORM_CX_MPI
   integer :: i, M0   

   M0 = M + 1
   VJ(1:N, M0) = WJ(1:N)   
  
   do I = 1, M
      dotwv = DOTPROD_SVCX_MPI(MPI_COMM_WORLD, N, VJ(1:N, I),VJ(1:N, M0), IERR)                
      VJ(1:N, M0) = VJ(1:N, M0) - dotwv * VJ(1:N, I)
   end do
  
   normw = NORM_CX_MPI(MPI_COMM_WORLD, N, VJ(1:N, M0), IERR)
   beta  = normw 

   if (normw == 0.0D0) then      
      LAN_MGS_ORTH_CX_MPI = .false.
   else
      VJ(1:N, M0) = VJ(1:N, M0) / normw
      LAN_MGS_ORTH_CX_MPI = .true.
   end if

end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             Grant-Schimit Orthogonization                  c
!c                  For the whole matrix                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

logical function GS_FULL_ORTH_CX_MPI(nRow, nCol, Mat, IERR)
   implicit none
   include 'mpif.h'
   integer, intent(IN)  :: nRow, nCol
   double complex, intent(INOUT)  :: Mat(nRow, nCol)
   integer, intent(out) :: ierr                         


!ccccccccccccccccc
   double complex     :: tmpV(nRow) 
   double precision   :: normV, NORM_CX_MPI
   logical            :: GS_ORTH_CX_MPI
   
   integer :: i   

   normV = NORM_CX_MPI(MPI_COMM_WORLD, nRow, Mat(1:nRow, 1), IERR)

   GS_FULL_ORTH_CX_MPI = .false.
   if (normV == 0.0D0)      return

   Mat(1:nRow,1) = Mat(1:nRow, 1)/normV  

   do I = 2, nCol      
       tmpV(1:nRow)  = Mat(1:nRow, I)
       if ( .not. GS_ORTH_CX_MPI(nRow, tmpV, I-1, Mat, ierr) )   return
   end do

   GS_FULL_ORTH_CX_MPI = .true.

end 
!****************************************************************


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Modified Grant-Schimit Orthogonization               c
!c                  For the whole matrix                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

logical function MGS_FULL_ORTH_CX_MPI(nRow, nCol, Mat, ierr)
   implicit none
   include 'mpif.h'
   integer, intent(IN)  :: nRow, nCol
   double complex, intent(INOUT)  :: Mat(nRow, nCol)
   integer, intent(out) :: ierr


!ccccccccccccccccc
   double complex   :: tmpV(nRow)
   double precision :: normV, NORM_CX_MPI  
   logical          :: MGS_ORTH_CX_MPI 

   integer :: i

   normV = NORM_CX_MPI(MPI_COMM_WORLD,nRow,Mat(1:nRow,1),IERR)

   MGS_FULL_ORTH_CX_MPI = .false.
   if (normV == 0.0D0)    return

   Mat(1:nRow,1) = Mat(1:nRow, 1)/normV

   do I = 2, nCol
       tmpV(1:nRow)  = Mat(1:nRow, I)
       if ( .not. MGS_ORTH_CX_MPI(nRow, tmpV, I-1, Mat, ierr) )  return
   end do

   MGS_FULL_ORTH_CX_MPI = .true.
   return

end 
!****************************************************************


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    Test whether a matrix is united orthogonal
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function isOrth_CX_MPI (nRow, nCol, A, epsi)
    implicit none
    include 'mpif.h'
    integer, intent(IN) :: nRow, nCol
    double complex, intent(IN)   :: A(nRow, nCol)
    double precision, intent(IN) :: epsi    

    double complex :: tmp

    integer :: I, J, ierr
    double complex, dimension(nCol, nCol) :: ATA

    isOrth_CX_MPI = .false.    

    call DOTPROD_MVCX_MPI(MPI_COMM_WORLD,nRow,nCol,A,nCol,A,ATA,ierr)
  
    do I = 1, nCol
       do J = 1, nCol         
           tmp = ATA(I,J)
           if (I == J) then          
               if (abs(tmp-1.0D0) > EPSI)    return
           else
               if (abs(tmp) > EPSI)          return
           end if
       end do
    end do           

end 
!cccccccccccccccccccccccccccccccccccccccc


