!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Create H Matrix in PIST:  H(I, J) = <VI|H|VJ>     c
!c    VI is the PIST vectors:                         c          
!c        WJ = (H-E)^-1*VI(I)                         c
!c        VI(I+1) = WJ - Sum((VI(K):WJ)*VI(K))        c
!c                 k = 1, 2, ..., I                   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function pistHij_SX_MPI(MYID,ROOTID,N,X,M,H0X_MPI,LinSolv_MPI,HMAT)
   implicit none
   include 'mpif.h'
   integer, intent(IN)   :: MYID, ROOTID
   integer, intent(IN)   :: N, M
   double complex, intent(IN)    :: X(N)
   external              :: H0X_MPI
   integer, external     :: LinSolv_MPI
   double complex,intent(OUT) :: HMAT(M,M)


   double complex :: VJ(N,M),WJ(N),TMPDOT(M)
   double complex :: TMP, DOT_CX_MPI, dot_cx
   double precision :: RES
   integer :: I, J, NUM, IERR
   logical :: GS_ORTH_SX_MPI
  
   tmp = DOT_CX_MPI(MPI_COMM_WORLD,N,X(1:N),X(1:N),IERR)
   tmp = sqrt(tmp)

   pistHIJ_SX_MPI = 0

   if (tmp == 0.0D0)        return

   VJ(1:N, 1)   = X(1:N)/tmp  

   PISTHIJ_SX_MPI = M
   do J = 1, M-1          ! WJ = VJ(1:N, J) / (H-EI)

     NUM = LinSolv_MPI(N, VJ(1:N, J), WJ, RES)
                                       
    if (NUM <= 0) then       
         PISTHIJ_SX_MPI = num
         return
    end if       

    if (.not. GS_ORTH_SX_MPI(N, WJ, J, VJ, IERR)) then        
         PISTHIJ_SX_MPI = -J;          return        
    end if    ! VJ(1:N, J+1) = H*VJ

  end do

  do J = 1, M
     call H0X_MPI(N, VJ(1:N, J), WJ, IERR)         ! WJ = H * VJ
  
     do I = 1, M
        TMPDOT(I)  = dot_cx(N, VJ(1:N,I), WJ(1:N))  ! HMAT(I,J) = [UI|H|UJ]
     end do

     call MPI_REDUCE(TMPDOT(1:M), HMAT(1:M, J), M, MPI_DOUBLE_COMPLEX,  &
                         MPI_SUM, ROOTID, MPI_COMM_WORLD, IERR)
  end do

end 


!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Create H Matrix in PIST:  H(I, J) = <VI|H|VJ>      c
!c    VI is the PIST vectors:                         c          
!c        WJ = (H-E)^-1*VI(I)                         c
!c        VI(I+1) = WJ - Sum((VI(K):WJ)*VI(K))        c
!c                 k = 1, 2, ..., I                   c
!c   This version need to store the old vectors and   c
!c   HMAT, but it will save the time to calculate the c
!c   previous calculated HIJ and vectors.             c
!c   The previous calculated ones are HMAT(1:M_LOW-1, c
!      1:M_LOW-1), and VJ(N, 1:M_LOW-1)               C
!c                                                    c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function pistH0_SX_MPI(MYID, ROOTID, N, X, M_MAX,  &
             M_LOW, M_HIGH, H0X_MPI, LinSolv_MPI,  VJ, HMAT)
   implicit none
   include 'mpif.h'
   integer, intent(IN)   :: MYID, ROOTID
   integer, intent(IN)   :: N, M_MAX, M_LOW, M_HIGH
   double complex, intent(IN)    :: X(N)
   external             :: H0X_MPI
   integer, external    :: LinSolv_MPI
   double complex,intent(INOUT) :: VJ(N, M_HIGH)
   double complex,intent(INOUT) :: HMAT(M_MAX, M_HIGH)


   double complex :: WJ(N), TMPDOT(M_HIGH)
   double complex :: TMP, DOT_CX_MPI,dot_cx
   double precision :: RES
   integer :: I, J, NUM, M_START, IERR
   logical :: GS_ORTH_SX_MPI
  
   pistH0_SX_MPI= 0

   if ((M_LOW<1) .or. (M_LOW>=M_HIGH) .or. (M_HIGH>M_MAX) )  return

   if (M_LOW == 1) then
      M_START = 1

      tmp = DOT_CX_MPI(MPI_COMM_WORLD,N,X(1:N),X(1:N),IERR)
      tmp = sqrt(tmp)

      if (tmp == 0.0D0)        return

      VJ(1:N, 1)   = X(1:N)/tmp  
   else
      M_START = M_LOW - 1
   end if

   PISTH0_SX_MPI = M_HIGH
   do J = M_START, M_HIGH-1          ! WJ = VJ(1:N, J) / (H-EI)
     NUM = LinSolv_MPI( N, VJ(1:N, J), WJ, RES)
                                       
    if (NUM <= 0) then       
         PISTH0_SX_MPI = num;  return
    end if       

    if (.not. GS_ORTH_SX_MPI(N, WJ, J, VJ, IERR)) then        
         PISTH0_SX_MPI = -J;        return        
    end if    ! VJ(1:N, J+1) = H*VJ

  end do
 
  do J = M_LOW, M_HIGH

     call H0X_MPI(N, VJ(1:N, J), WJ, IERR)         ! WJ = H * VJ

     do I = 1, M_HIGH
        TMPDOT(I) = dot_cx(N, VJ(1:N,I), WJ(1:N)) 
     end do

     call MPI_REDUCE(TMPDOT(1:M_HIGH), HMAT(1:M_HIGH, J), M_HIGH, MPI_DOUBLE_COMPLEX,   &
                         MPI_SUM, ROOTID, MPI_COMM_WORLD, IERR)
  end do

  do J = 1, M_LOW-1

     call H0X_MPI(N, VJ(1:N, J), WJ, IERR)         ! WJ = H * VJ

     do I = M_LOW, M_HIGH
        TMPDOT(I) = dot_CX(N, VJ(1:N,I), WJ(1:N)) 
     end do

     call MPI_REDUCE(TMPDOT(M_LOW:M_HIGH), HMAT(M_LOW:M_HIGH, J), (M_HIGH-M_LOW+1),     &
                  MPI_DOUBLE_COMPLEX,MPI_SUM, ROOTID, MPI_COMM_WORLD, IERR)
  end do

end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




