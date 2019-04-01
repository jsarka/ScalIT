!ccccccccccccccccccccccccccccccccccccccc
!c  
!c        testing program for PIST
!c
!ccccccccccccccccccccccccccccccccccccccc

      program lan_dx_mpi_drv

      INCLUDE 'mpif.h'
      integer,parameter  :: N = 100
      INTEGER, PARAMETER :: N1= 25
      integer, parameter :: EIGNUM = 4
      double precision, parameter :: E0     = -10.0D0
      double precision, parameter :: ETOL   = 1.0e-7
      double precision :: ct1, ct2

      double precision, dimension(EIGNUM) ::  EIG
      double complex, dimension(N1)      ::  x   
      external :: hmult, pmult     

      integer, parameter :: ROOTID = 0
      integer  :: nodnum, myid,ierr
      integer  ::  i,iter, qmr, num  
      double precision :: res, MT1, MT2
      integer  :: testType ;    ! testing
      integer  :: MAX_NUM  = 100

      CALL MPI_INIT(IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NODNUM, IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)
     
      call HINIT_CX()
!      res = max(X(1:N))


      if (myid==rootid) then
        print *, ' test for Lan subroutine'
        print *, ' Input testing type: '
        print *, ' -1/0: for DSYEV directly'
        print *, ' 1,2:  for fast LanCONV '
        print *, ' 3,4:  for slow LanCONV '     
        read(*, *) testType
      end if

      call MPI_BCAST(testType, 1, mpi_integer, rootid, mpi_comm_world, ierr)

      NUM = MYID*N1
      do i=1,n1
         x(i) = (i+NUM) * 1.0d0 
      end do

     call cpu_time(ct1)
     MT1 = MPI_WTIME()

     IF (MYID==ROOTID) THEN
         PRINT *, '# OF NODE: ', NODNUM
     END IF

!     if(myid == 1)  print *, 'X for ID=', myID, x
 
    select case (testType)
    case (-1)
      IF (MYID==ROOTID) THEN
          print *
          print *, '  Get Part of Eigen values directly'
          print *, ' Input number of eigen values:'      
          read *, num
          call printEig0_CX(E0, eignum)
      END IF

    case (0)
      IF (MYID==ROOTID) THEN
          print *
          print *, 'Get Eigen values directly:'      
          call printEig_CX()
      END IF

     case(1)      
      if(lanConv1_DX_MPI(MYID, ROOTID, E0, ETOL, N1, X, EIGNUM, MAX_NUM,  &
                       HMULT, EIG, RES)>0) THEN
         IF (MYID==ROOTID) THEN
           print *
           print *, ' Performing Lanczos using faster Convergent Criterion'
           print *, 'Res:', res
           print *, 'Eig: ', eig
         END IF
      else
        IF (MYID==ROOTID) THEN
           print *, 'Error in LanConv1_MPI'
           print *, 'Res:', res
           print *, 'Eig: ', eig
        END IF
      end if

     case(2)
      if(lanConv2_dx_MPI(MYID, ROOTID,E0, ETOL, N1, X, EIGNUM, MAX_NUM,  &
                       HMULT, EIG, RES)>0) THEN
          IF (MYID==ROOTID) THEN
              print *
              print *, ' Performing Lanczos using Faster Convergent Criterion'
              print *, 'Res:', res
              print *, 'Eig: ', eig
          END IF
      else
         IF (MYID==ROOTID) THEN
             print *, 'Error in LanConv2_MPI'
             print *, 'Res:', res
             print *, 'Eig: ', eig
          END IF
      end if

     case(3)
      if(lanConv51_dx_MPI(MYID, ROOTID,E0, ETOL, N1, X, EIGNUM, MAX_NUM,  &
                       HMULT, EIG, RES)>0) THEN
         IF (MYID==ROOTID) THEN
              print *
              print *, ' Performing less strict Lanczos using Slower Convergent Criterion'
              print *, 'Res:', res
              print *, 'Eig: ', eig
          END IF
      else
         IF (MYID==ROOTID) THEN
            print *, 'Error in LanConv51_MPI'
            print *, 'Res:', res
            print *, 'Eig: ', eig
         END IF
      end if

    case (4)
      if(lanConv52_dx_MPI(MYID, ROOTID, E0, ETOL, N1, X, EIGNUM, MAX_NUM,  &
                       HMULT, EIG, RES)>0) THEN
         IF (MYID==ROOTID) THEN
            print *
            print *, ' Performing less strict Lanszos using Slower Convergent Criterion'
            print *, 'Res:', res
            print *, 'Eig: ', eig
         END IF
      else
        IF (MYID==ROOTID) THEN
          print *, 'Error in LanConv52_MPI'
          print *, 'Res:', res
          print *, 'Eig: ', eig
        END IF
      end if

    end select 
    call cpu_TIME(ct2)
    MT2 = MPI_WTIME()

    IF (MYID == ROOTID) THEN
      print *, 'CPU Time:', ct2-ct1
      PRINT *, 'MPI TIME:', MT2-MT1
      print *, 'Finish Lanczos testing'
    END IF
    CALL MPI_FINALIZE(IERR)
     
    END

!c
!c     Get symmetric h matrix
!c
      SUBROUTINE  HINIT_cx()
      
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
            H(I, J)  = DCMPLX(rand(), rand())
            H(J, I)  = DCONJG(H(I, J))
         end do
      end do
      do I=1,N
          H(I,I) = rand()
      end do

      END IF

      CALL MPI_BCAST(H, N*N, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, IERR)
      end

!c
!c     y = H x
!c
      SUBROUTINE  HMULT(N1, X0, Y, ierr)
      IMPLICIT NONE
      integer,parameter :: N = 100
      integer :: N1  
      double complex, dimension(N1)   :: X0
      DOUBLE complex, DIMENSION(N1)   :: Y
      double complex, dimension(n, n) ::  H

      integer :: MYID, IERR, SIZE
      INCLUDE 'mpif.h'
      double complex , dimension(N) :: X

      common /H/ H
    
      integer :: i, j,K, INDEX, M

      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, SIZE, IERR)
      call MPI_ALLGATHER(x0,n1,MPI_DOUBLE_COMPLEX, X, N1, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr )

!      IF (MYID == 0) THEN
!          print *, 'call HMult, x0:',X0
!          print *, 'All X:', X
!      END IF

      INDEX = N/SIZE * MYID

      do i = 1, N1
         y(i) = 0.0d0
         K = I+ INDEX
         do j=1, N
            y(i) = y(i) + H(K,j)*x(j)
         end do
      end do 
        
!      IF (MYID == 0) THEN
!          print *, 'Finish HMult', y
!      END IF
 
      end

!c
!c     y = P^(-1) x, Jacobi preconditioning
!c
      SUBROUTINE  PMULT(m, X, Y, ierr)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      integer,parameter :: N = 100
      INTEGER, PARAMETER :: N1 = 25

      double precision, dimension(n)    :: X
      DOUBLE PRECISION, DIMENSION(N1)   :: Y
      double precision, dimension(n, n) ::  H

      common /H/  H
    
      integer :: i, MYID, SIZE, INDEX, M, IERR,J


      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, SIZE, IERR)
   
!      IF (MYID == 0) THEN
!          print *, 'call PMult'
!      END IF

      INDEX = N/SIZE * MYID
      do  i = 1, N1
           J = INDEX + I
         y(I) = x(J)/H(J,J)
      end do
         
      end



!c
!c   Get the eigen values of H
!c
      subroutine printEig0_cx(E0, M)
      integer,parameter :: N = 100

      double precision, intent(IN)  :: E0
      integer, intent(IN)           :: M

      double precision, dimension(n)    :: allEig
      double complex, dimension(n, n)   ::  H
      double complex, dimension(3*n)    :: work
      double precision, dimension(M)    :: eig
      double precision, dimension(3*n)  :: rwork
      integer :: lwork = 3*N
      integer :: tmp
  
      common /H/  H
      
      CALL ZHEEV('N', 'U', N, H, N, allEig, work, lwork, rwork, tmp)
      
      if (tmp == 0) then
          print *, M,' Eigen values of H around ', E0
          call getWindow(E0, N, allEig, M, eig)
          call Reorder('A', M, eig)
          print *, eig
      else
          print *, 'Error while solving eigen value problem of H'
      end if

      end

      subroutine printEig_CX()
      integer,parameter :: N = 100

      double precision, dimension(n)    :: allEig
      double complex, dimension(n, n)   ::  H
      double precision, dimension(3*n)  :: work, rwork
      double precision, dimension(N)    :: eig
      integer :: lwork = 3*N
      integer :: tmp     
  
      common /H/  H
      

!      print *, 'Matrix'
!      print *, H

      CALL ZHEEV('N', 'U', N, H, N, allEig, work, lwork, rwork, tmp)
      
      if (tmp == 0) then
          print *, ' Eigen values of H around ', allEig
      else
          print *, 'Error while solving eigen value problem of H'
      end if

      end

