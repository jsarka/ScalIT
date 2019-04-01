!ccccccccccccccccccccccccccccccccccccccc
!c  
!c        testing program for QMR
!c
!ccccccccccccccccccccccccccccccccccccccc

      program qmr_mpi_drv
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      integer,parameter  :: N = 100
      INTEGER, PARAMETER :: N1 = 25
      integer, parameter :: QMRMAX  = 1000
      INTEGER, PARAMETER :: QMRTYPE = 0
      double precision, parameter :: QMRTOL = 1.0D-4

      double precision :: ct1, ct2

      double precision, dimension(N)      ::  x0   
      DOUBLE PRECISION, DIMENSION(N1)     ::  X, B
      external :: hmult, pmult     

      integer, parameter :: ROOTID = 0
      integer  :: nodnum, myid,ierr
      integer  ::  i,iter, qmr_MPI, num  
      double precision :: res, MT1, MT2

      CALL MPI_INIT(IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NODNUM, IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)
     
      call HINIT()
!      res = max(X(1:N))

      NUM = MYID*N1
      DO I = 1, N1 
         B(I) = (I+NUM)*1.0D0
      END DO
!      do i=1,n
!         x0(i) = i * 1.0d0 
!      end do

     call cpu_time(ct1)
     MT1 = MPI_WTIME()

     IF (MYID==ROOTID) THEN
         print *, ' TEST QMR_MPI'
         PRINT *, '# OF NODE: ', NODNUM
         PRINT *, 'MY ID:', MYID
         PRINT *, 'ORIGINAL X0:'
         PRINT *, B 
     END IF
 
      NUM = QMR_MPI(QMRTYPE, QMRMAX, QMRTOL, N1, B,        &
              HMULT, PMULT, X, RES, IERR)
      call mpi_GATHER(X, N1, MPI_DOUBLE_PRECISION, X0, N1,            &
                MPI_DOUBLE_PRECISION, ROOTID, MPI_COMM_WORLD, IERR)
      call cpu_TIME(ct2)
      MT2 = MPI_WTIME()

    IF (MYID == ROOTID) THEN
      PRINT *, ' RESULT AT ROOT=', MYID, 'X:'
      PRINT *, X0
      print *, 'CPU Time:', ct2-ct1
      PRINT *, 'MPI TIME:', MT2-MT1
      print *, 'Finish PIST testing'
    END IF
    CALL MPI_FINALIZE(IERR)
     
      END

!c
!c     Get symmetric h matrix
!c
      SUBROUTINE  HINIT()
      
      integer,parameter :: N  = 100 
      double precision, dimension(N, N):: H                    

      common /H/  H

      integer :: i, j
      integer :: MYID, IERR
      INCLUDE 'mpif.h'
      
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)
      if (MYID == 0) THEN
      do I=1,N
         do J = 1, I
            H(I, J)  = rand()
            H(J, I)  = H(I, J)
         end do
      end do
!       PRINT *, 'CALL HINIT'    
      END IF

      CALL MPI_BCAST(H, N*N, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)
      end

!c
!c     y = H x
!c
      SUBROUTINE  HMULT(m, X, Y, ierr)
      IMPLICIT NONE
      integer,parameter :: N = 100
      integer,parameter :: N1 = 25
      INTEGER :: M
      double precision, dimension(n1)    :: X
      DOUBLE PRECISION, DIMENSION(N1)   :: Y
      double precision, dimension(n, n) ::  H

      integer :: MYID, IERR, SIZE
      INCLUDE 'mpif.h'
      double precision, dimension(N) :: xtmp

      common /H/ H
    
      integer :: i, j,K, INDEX

      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, SIZE, IERR)

      call MPI_ALLGATHER(X, N1, MPI_DOUBLE_precision, XTMP, N1, &
                  MPI_DOUBLE_precision, MPI_COMM_WORLD, IERR )
   

!      IF (MYID == 0) THEN
!          print *, 'call HMult'
!      END IF

      INDEX = N/SIZE * MYID

      do i = 1, N1
         y(i) = 0.0d0
         K = I+ INDEX
         do j=1, N
            y(i) = y(i) + H(K,j)*xtmp(j)
         end do
      end do 
         
      end

!c
!c     y = P^(-1) x, Jacobi preconditioning
!c
      SUBROUTINE  PMULT(m, X, Y, ierr)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      integer,parameter :: N = 100
      INTEGER, PARAMETER :: N1 = 25

      INTEGER :: M
      double precision, dimension(n1)    :: X
      DOUBLE PRECISION, DIMENSION(N1)   :: Y
      double precision, dimension(n, n) ::  H

      common /H/  H
    
      integer :: i,J, MYID, SIZE0, INDEX, IERR
      double precision, dimension(N) ::xtmp

      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, SIZE0, IERR)

      call MPI_ALLGATHER(X, N1, MPI_DOUBLE_precision, XTMP, N1, &
                  MPI_DOUBLE_precision, MPI_COMM_WORLD, IERR )   
   
!      IF (MYID == 0) THEN
!          print *, 'call PMult'
!      END IF

      INDEX = N/SIZE0 * MYID
      do  i = 1, N1
           J = INDEX + I
         y(I) = xtmp(J)/H(J,J)
      end do
         

      end
