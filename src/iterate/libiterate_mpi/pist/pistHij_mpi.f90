!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Create H Matrix in PIST: H(I, J) = <VI|H|VJ>      c
!c    VI is the PIST vectors:                         c          
!c        WJ = (H-E)^-1*VI(I)                         c
!c        VI(I+1) = WJ - Sum((VI(K):WJ)*VI(K))        c
!c                 k = 1, 2, ..., I                   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function pistHij_MPI(MYID,ROOTID,N,X,M,H0X_MPI,LinSolv_MPI,HMAT)
   implicit none
   include 'mpif.h'
   integer, intent(IN)   :: MYID, ROOTID
   integer, intent(IN)   :: N, M
   double precision, intent(IN)    :: X(N)
   external              :: H0X_MPI
   integer, external     :: LinSolv_MPI
   double precision,intent(OUT) :: HMAT(M,M)


   double precision :: VJ(N,M),WJ(N),TMPDOT(M),SUMDOT(M)
   double precision :: TMP, RES, NORM_MPI
   integer :: I, J, NUM, IERR
   logical :: MGS_ORTH_MPI
  
   tmp  = NORM_MPI(MPI_COMM_WORLD, N, X(1:N), IERR)

   pistHIJ_MPI = 0

   if (tmp == 0.0D0)    return

   VJ(1:N, 1)   = X(1:N)/tmp  

   PISTHIJ_MPI = M
   do J = 1, M-1          ! WJ = VJ(1:N, J) / (H-EI)

     NUM = LinSolv_MPI(N, VJ(1:N, J), WJ, RES)
                                       
    if (NUM <= 0) then       
         PISTHIJ_MPI = num;      return
    end if       

    !if (.not. GS_ORTH_MPI(N, WJ, J, VJ, IERR)) then        
    if (.not. MGS_ORTH_MPI(N, WJ, J, VJ, IERR)) then        
         PISTHIJ_MPI = -J;       return        
    end if    ! VJ(1:N, J+1) = H*VJ

  end do

  do J = 1, M

     call H0X_MPI(N, VJ(1:N, J), WJ, IERR)         ! WJ = H * VJ
  
     do I = J, M
        TMPDOT(I)  = dot_product(VJ(1:N,I), WJ(1:N))  ! HMAT(I,J) = <UI|H|UJ>
     end do

     call MPI_REDUCE(TMPDOT(J:M), HMAT(J:M, J), (M-J+1), MPI_DOUBLE_PRECISION,  &
                         MPI_SUM, ROOTID, MPI_COMM_WORLD, IERR)

     if (MYID == ROOTID) then
        do I = J, M       
           HMAT(J, I)  = HMAT(I, J)
        end do
     end if
  end do

end 


!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Create H Matrix in PIST,  H(I, J) = <VI|H|VJ>     c
!c    VI is the PIST vectors:                         c          
!c        WJ = (H-E)^-1*VI(I)                         c
!c        VI(I+1) = WJ - Sum((VI(K):WJ)*VI(K))        c
!c                 k = 1, 2, ..., I                   c
!c   This version need to store the old vectors and   c
!c   HMAT, but it will save the time to calculate the c
!c   previous calculated HIJ and vectors.             c
!c   The previous calculated ones are HMAT(1:M_LOW-1, c
!      1:M_LOW-1), and VJ(N, 1:M_LOW-1)               C
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function pistH0_MPI(MYID, ROOTID, N, X, M_MAX,   &
             M_LOW, M_HIGH, H0X_MPI, LinSolv_MPI, VJ, HMAT, num2)
   implicit none
   include 'mpif.h'
   integer, intent(IN)   :: MYID, ROOTID
   integer, intent(IN)   :: N, M_MAX, M_LOW, M_HIGH
   double precision, intent(IN)    :: X(N)
   external             :: H0X_MPI
   integer, external    :: LinSolv_MPI
   double precision, intent(INOUT) :: VJ(N, M_HIGH)
   double precision, intent(INOUT) :: HMAT(M_MAX, M_HIGH)


   double precision :: WJ(N)
   double precision :: TMPDOT(M_HIGH)
   double precision :: TMP, RES, NORM_MPI
   integer :: I, J, NUM, M_START, IERR, num2
   logical :: MGS_ORTH_MPI
 
   double precision :: wtp1, wtp2, wtp3, wtp4
   double precision :: wtp5, wtp6, wtp7, wtp8
   double precision :: wtp9, wtp10, wtp11, wtp12
 
   pistH0_MPI= 0

   if ((M_LOW<1) .or. (M_LOW>=M_HIGH) .or. (M_HIGH>M_MAX) )  return

   if (M_LOW == 1) then
      M_START = 1

      tmp  = NORM_MPI(MPI_COMM_WORLD, N, X(1:N), IERR)

      if (tmp == 0.0D0)         return

      VJ(1:N, 1)   = X(1:N)/tmp  
   else
      M_START = M_LOW - 1
   end if

   num2 = 0
   PISTH0_MPI = M_HIGH
   do J = M_START, M_HIGH-1          ! WJ = VJ(1:N, J) / (H-EI)

     wtp1 = MPI_WTime()

     NUM = LinSolv_MPI(N, VJ(1:N, J), WJ, RES)

     if (myid==rootid)              &
        print *,'# of Lanczos iter.:',J,'. # of iter. for Linear Solver:',num

     num2 = num2 + num

     wtp2 = MPI_WTime()
     if (myid==rootID) print *, ' MPI Time for QMR:',wtp2-wtp1
     if (myid==rootID) print *
                                      
    if (NUM <= 0) then       
         PISTH0_MPI = num;       return
    end if       

    !if (.not. GS_ORTH_MPI(N, WJ, J, VJ, IERR)) then        
    if (.not. MGS_ORTH_MPI(N, WJ, J, VJ, IERR)) then        
         PISTH0_MPI = -J;        return        
    end if    ! VJ(1:N, J+1) = H*VJ

    !wtp3 = MPI_WTime()
    !if (myid==rootID) print *, ' MPI Time for ORTH:',wtp3-wtp2
    !if (myid==rootID) print *

   end do

   wtp4 = MPI_WTime() 

  do J = M_LOW, M_HIGH

     !wtp6 = MPI_WTime()

     call H0X_MPI(N, VJ(1:N, J), WJ, IERR)         ! WJ = H * VJ

     !wtp7 = MPI_WTime()
     !if (myid==rootID) print *, ' MPI Time for H0X:',wtp7-wtp6

     do I = 1, M_HIGH
        TMPDOT(I) = dot_product(VJ(1:N,I), WJ(1:N)) 
     end do

     !wtp8 = MPI_WTime()
     !if (myid==rootID) print *, ' MPI Time for DP :',wtp8-wtp7

     call MPI_REDUCE(TMPDOT(1:M_HIGH), HMAT(1:M_HIGH, J), M_HIGH, MPI_DOUBLE_PRECISION,   &
                         MPI_SUM, ROOTID, MPI_COMM_WORLD, IERR)
     !wtp9 = MPI_WTime()
     !if (myid==rootID) print *, ' MPI Time for Red:',wtp9-wtp8
     !if (myid==rootID) print *

     if (MYID == ROOTID) then
        do I = 1, M_HIGH    
           HMAT(J, I)  = HMAT(I, J)
        end do
     end if

  end do

  wtp5 = MPI_WTime()
  if (myid==rootID) print *, ' MPI Time for HMAT:',wtp5-wtp4
  !if (myid==rootID) print *

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




