!ccccccccccccccccccccccccccccccccccccccc
!c  
!c        testing program for PIST
!c
!ccccccccccccccccccccccccccccccccccccccc

      program pist_cx_mpi_drv
      implicit none
      INCLUDE 'mpif.h'
      integer,parameter  :: N = 100
      INTEGER, PARAMETER :: N1= 25
      integer, parameter :: EIGNUM = 4
      double precision, parameter :: E0     = 0.0D0
      double precision, parameter :: ETOL   = 1.0e-4
      double precision :: ct1, ct2

      double complex, dimension(EIGNUM) ::  EIG
      double complex, dimension(N1)      ::  x   
      external :: hmult
      integer, external :: mySolvCX     
      integer  :: pist1_cx_mpi, pist2_cx_mpi,pist51_cx_mpi,pist52_cx_mpi
      logical  :: pistconverg

      integer, parameter :: ROOTID = 0
      integer  :: nodnum, myid,ierr
      integer  ::  i, num  
      double precision :: res, MT1, MT2
      integer  :: testType = 0;
      integer  :: MAX_NUM  = 50
 
      CALL MPI_INIT(IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NODNUM, IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)
     
      call HINIT()
!      res = max(X(1:N))

      IF (MYID == ROOTID) THEN
         print *, ' Test for PIST subroutine'
         print *, ' Input testing type: '
         print *, '-1/0: for ZGEEV directly'
         print *, ' 1,2:  for fast PISTCONV '
         print *, ' 3,4:  for slow PISTCONV '     
         read(*, *) testType 
      END IF

      CALL MPI_BCAST(TESTTYPE, 1, MPI_INTEGER, ROOTID, MPI_COMM_WORLD, IERR)

!     TESTTYPE = 1
      NUM = MYID*N1
      do i=1,n1
         x(i) = (i+NUM) * 1.0d0 
      end do

     call cpu_time(ct1)
     MT1 = MPI_WTIME()

     IF (MYID==ROOTID) THEN
         PRINT *, '# OF NODE: ', NODNUM
         PRINT *, 'root ID:', MYID
     END IF
 
     select case (testType)
     case (-1)
       IF (MYID==ROOTID) THEN
          print *
          print *, '  Get Part of Eigen values directly'
          print *, 'Input the number of eigen values:'
          read(*, *) num
          call printEig0_CX(E0, num)
       END IF

     case (0)
       IF (MYID==ROOTID) THEN
          print *
          print *, 'Get Eigen values directly:'      
          call printEig_CX()
       END IF

     case(1)      
      if(pist1_cx_MPI(MYID, ROOTID, E0, ETOL, N1, X, EIGNUM, MAX_NUM,  &
                       HMULT,  MySolvCX, EIG, RES)>0) THEN
         IF (MYID==ROOTID) THEN
           print *
           print *, ' Performing PIST using faster Convergent Criterion'
           print *, 'Res:', res
           print *, 'Eig: ', eig
         END IF
      else
        IF (MYID==ROOTID) THEN
           print *, 'Error in PISTCONV'
        END IF
      end if

     case(2)
      if(pist2_CX_MPI(MYID, ROOTID,E0, ETOL, N1, X, EIGNUM, MAX_NUM,  &
                       HMULT,  MySolvCX, EIG, RES)>0) THEN
          IF (MYID==ROOTID) THEN
              print *
              print *, ' Performing PIST using Faster Convergent Criterion'
              print *, 'Res:', res
              print *, 'Eig: ', eig
          END IF
      else
         IF (MYID==ROOTID) THEN
             print *, 'Error in PISTCONVERG'
          END IF
      end if

     case(3)
      if(pist51_CX_MPI(MYID, ROOTID,E0, ETOL, N1, X, EIGNUM, MAX_NUM,  &
                       HMULT,  MySolvCX, EIG, RES)>0) THEN
         IF (MYID==ROOTID) THEN
              print *
              print *, ' Performing less strict PIST using Slower Convergent Criterion'
              print *, 'Res:', res
              print *, 'Eig: ', eig
          END IF
      else
         IF (MYID==ROOTID) THEN
            print *, 'Error in PISTCONVERG'
         END IF
      end if

    case (4:)
      if(pist52_CX_MPI(MYID, ROOTID, E0, ETOL, N1, X, EIGNUM, MAX_NUM,  &
                       HMULT,  MySolvCX, EIG, RES)>0) THEN
         IF (MYID==ROOTID) THEN
            print *
            print *, ' Performing less strict PIST using Slower Convergent Criterion'
            print *, 'Res:', res
            print *, 'Eig: ', eig
         END IF
      else
        IF (MYID==ROOTID) THEN
          print *, 'Error in PISTCONVERG'
        END IF
      end if

    end select 
    call cpu_TIME(ct2)
    MT2 = MPI_WTIME()

    IF (MYID == ROOTID) THEN
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
      double complex, dimension(N, N):: H                    

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
       PRINT *, 'CALL HINIT'    

      do I = 1, N
          H(I, I) = DCMPLX(rand(), rand())
      end do

      END IF

      CALL MPI_BCAST(H, N*N, MPI_DOUBLE_complex, 0, MPI_COMM_WORLD, IERR)
      end

!c
!c     y = H x
!c
      SUBROUTINE  HMULT(N1, X0, Y, ierr)
      IMPLICIT NONE
      integer,parameter :: N = 100
      integer,intent(IN) :: N1 
      double COMPLEX, dimension(n1)   :: X0
      DOUBLE COMPLEX, DIMENSION(N1)   :: Y
      double COMPLEX, dimension(n, n) ::  H

      integer :: MYID, IERR, SIZE
      INCLUDE 'mpif.h'
      double complex, dimension(N)  :: X      

      common /H/ H
    
      integer :: i, j,K, INDEX, M

      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, SIZE, IERR)

      call MPI_ALLGATHER(x0,n1,MPI_DOUBLE_COMPLEX, X, N1, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr )

!      IF (MYID == 0) THEN
!          print *, 'call HMult'
!      END IF

      INDEX = N/SIZE * MYID

      do i = 1, N1
         y(i) = 0.0d0
         K = I+ INDEX
         do j=1, N
            y(i) = y(i) + H(K,j)*x(j)
         end do
      end do 
         
      end

!c
!c     y = P^(-1) x, Jacobi preconditioning
!c
      SUBROUTINE  PMULT(N1, X0, Y, ierr)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      integer,parameter :: N = 100
      INTEGER, intent(IN) :: N1 
      double COMPLEX, dimension(N1)   :: X0
      DOUBLE COMPLEX, DIMENSION(N1)   :: Y
      double COMPLEX, dimension(n, n) ::  H

      common /H/  H
    
      integer :: i, MYID, SIZE, INDEX, M, IERR,J
      double complex, dimension(N) :: X

      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, SIZE, IERR)  
      call MPI_ALLGATHER(x0,n1,MPI_DOUBLE_COMPLEX, X, N1, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr )
!      IF (MYID == 0) THEN
!          print *, 'call PMult'
!      END IF

      INDEX = N/SIZE * MYID
      do  i = 1, N1
           J = INDEX + I
         y(I) = x(J)/H(J,J)
      end do
         
      end



   integer function mySolvCX(N1, X, Y, res)
      implicit none
      integer, intent(IN) :: N1
      double complex, intent(IN), dimension(N1) :: X
      double complex, intent(OUT),dimension(N1) :: Y
      double precision, intent(OUT) :: res

      integer, parameter :: QMRTYPE=0
      integer, parameter :: QMRMAX = 1000
      double precision, parameter :: QMRTOL = 1.0D-4
      integer :: QMR_CX_MPI
      external :: HMULT, PMULT
      
      mySolvCX = QMR_CX_MPI(QMRTYPE,QMRMAX,QMRTOL,N1,X,HMULT,PMULT,Y,RES)

   end 


!c
!c   Get the eigen values of H directly using ZGEEV
!c
      subroutine printEig_CX()
      integer,parameter :: N = 100
      integer, parameter :: LDVL=1, LDVR=1
      double complex :: VL, VR
      double complex, dimension(n)    :: eig
      double complex, dimension(n, n) ::  H
      double complex, dimension(3*n)  :: work
      double precision, dimension(2*n)  :: RWORK
      integer :: lwork = 3*N
      integer :: tmp
  
      common /H/  H

!      print *, 'Matrix:'
!      print *, H

      ! CALL DSYEV('N', 'U', N, H, N, Eig, work, lwork, tmp)
      CALL ZGEEV('N', 'N', N, H, N, Eig, VL, LDVL, VR, LDVR,   &
                 WORK, LWORK, RWORK, tmp)

   
      if (tmp == 0) then
          print *, 'Eigen values of H using ZGEEV'
          call Reorder_CX('A', 1, N, eig)
          print *, eig
      else
          print *, 'Error while solving eigen value problem of H'
      end if

      end

!c
!c   Get the eigen values of H directly using ZGEEV
!c
      subroutine printEig0_CX(E0, num )
      integer,parameter :: N = 100
      double precision, intent(IN) :: E0
      integer, intent(IN)  :: num

      integer, parameter :: LDVL=1, LDVR=1
      double complex :: VL, VR
      double complex, dimension(n)    :: eig
      double complex, dimension(n, n) ::  H
      double complex, dimension(3*n)  :: work
      double precision, dimension(2*n)  :: RWORK
      integer :: lwork = 3*N
      integer :: tmp, i
      double complex, dimension(num) :: eig0

      common /H/  H

!      print *, 'Matrix:'
!      print *, H

      ! CALL DSYEV('N', 'U', N, H, N, Eig, work, lwork, tmp)
      CALL ZGEEV('N', 'N', N, H, N, Eig, VL, LDVL, VR, LDVR,   &
                 WORK, LWORK, RWORK, tmp)
  
      if (tmp == 0) then
          print *, 'Eigen values of H using ZGEEV around ', E0
          call getWindow_CX(E0, 1, N, eig, num, eig0)
          call Reorder_CX('A', 1, num, eig0)
          print *, eig0
      else
          print *, 'Error while solving eigen value problem of H'
      end if

      end
