!ccccccccccccccccccccccccccccccccccccccc
!c  
!c        testing program for PIST
!c
!ccccccccccccccccccccccccccccccccccccccc

      program pist_drv
      implicit none
      integer,parameter  :: N = 100
      integer, parameter :: EIGNUM = 4
      double precision, parameter :: E0     = 0.0D0
      double precision, parameter :: ETOL   = 1.0e-4
      double precision :: ct1, ct2

      double precision, dimension(EIGNUM) ::  EIG
      double precision, dimension(N)      ::  x   
      integer  :: pist1, pist2, pist51, pist52
      external :: HMULT, QMRMULT

      integer  ::  i, num  
      double precision :: res
      integer  :: testType = 0;
      integer  :: MAX_NUM  = 50

      print *, ' test for PIST subroutine'
      print *, ' Input testing type: '
      print *, ' -1/0: for DSYEV directly'
      print *, ' 1,2:  for fast PISTCONV '
      print *, ' 3,4:  for slow PISTCONV '     
      read(*, *) testType
 
      call HINIT()
!      res = max(X(1:N))

      do i=1,n
         x(i) = i * 1.0d0 
      end do

     call cpu_time(ct1)
     select case (testType)
    case (-1)
      print *
      print *, '  Get Part of Eigen values directly'
      print *, ' Input number of eigen values:'      
      read *, num
      call printEig0(E0, num)

    case (0)
      print *
      print *, 'Get Eigen values directly:'      
      call printEig()

     case(1)
     print *
      print *, ' Performing PIST using faster Convergent Criterion'
      if(pist1(E0, ETOL, N, X, EIGNUM, MAX_NUM, HMULT,QMRMULT, EIG, RES)>0) THEN
         print *, 'Res:', res
         print *, 'Eig: ', eig
      else
         print *, 'Error in PISTCONV'
      end if

     case(2)
     print *
      print *, ' Performing PIST using Faster Convergent Criterion'
      if(pist2(E0, ETOL, N, X, EIGNUM, MAX_NUM, HMULT,QMRMULT, EIG, RES)>0) THEN
         print *, 'Res:', res
         print *, 'Eig: ', eig
      else
         print *, 'Error in PISTCONVERG'
      end if

     case(3)
     print *
      print *, ' Performing less strict PIST using Slower Convergent Criterion'
      if(pist51(E0, ETOL, N, X,EIGNUM, MAX_NUM, HMULT, QMRMULT, EIG, RES)>0) THEN
         print *, 'Res:', res
         print *, 'Eig: ', eig
      else
         print *, 'Error in PISTCONVERG'
      end if


    case (4)
      print *
      print *, ' Performing less strict PIST using Slower Convergent Criterion'
      if(pist52(E0, ETOL, N, X, EIGNUM, MAX_NUM, HMULT, QMRMULT,EIG, RES)>0) THEN
         print *, 'Res:', res
         print *, 'Eig: ', eig
      else
         print *, 'Error in PISTCONVERG'
      end if

    end select 
    call cpu_TIME(ct2)

      print *, 'CPU Time:', ct2-ct1
      print *, 'Finish PIST testing'

      END

!c
!c     Get symmetric h matrix
!c
      SUBROUTINE  HINIT()
      
      integer,parameter :: N  = 100 
      double precision, dimension(N, N):: H                    

      common /H/  H

      integer i, j
      
      do I=1,N
         do J = 1, I
            H(I, J)  = rand()
            H(J, I)  = H(I, J)
         end do
      end do    

      end

!c
!c     y = H x
!c
      SUBROUTINE  HMULT(m, X,Y)
      
      integer,parameter :: N = 100

      double precision, dimension(n)    :: X, Y
      double precision, dimension(n, n) ::  H

      common /H/ H
    
      integer i, j

!c      print *, 'call HMult'

      do i = 1, N
         y(i) = 0.0d0
         do j=1, N
            y(i) = y(i) + H(i,j)*x(j)
         end do
      end do 
         
      end

!c
!c     y = P^(-1) x, Jacobi preconditioning
!c
      subroutine PMULT(m, X, Y)

      integer,parameter :: N = 100

      double precision, dimension(n)    :: X, Y
      double precision, dimension(n, n) ::  H

      common /H/  H
    
      integer i
      
!c      print *, 'Call PMult'

      do  i = 1, N
         y(i) = x(i)/H(i,i)
      end do
         
      end

      integer function QMRMULT(M, B, X, RES)
      implicit none
      integer, parameter :: QMRMAX   = 1000
      integer, parameter :: QMR_TYPE = 0
      integer, parameter :: N = 100
      double precision, parameter :: QMRTOL = 1.0D-4      

      integer, intent(IN)  :: M
      double precision, dimension(M), intent(IN)  :: B
      double precision, dimension(M), intent(OUT) :: X
      double precision, intent(OUT) :: RES

      external :: HMULT, PMULT
      integer, external :: QMR
      
      QMRMULT = QMR(QMR_TYPE, QMRMAX,QMRTOL, N, B, HMULT, PMULT, X, RES)

      end

!c
!c   Get the eigen values of H
!c
      subroutine printEig0(E0, M)
      integer,parameter :: N = 100

      double precision, intent(IN)  :: E0
      integer, intent(IN)           :: M

      double precision, dimension(n)    :: allEig
      double precision, dimension(n, n) ::  H
      double precision, dimension(3*n)  :: work
      double precision, dimension(M)    :: eig
      integer :: lwork = 3*N
      integer :: tmp
  
      common /H/  H
      
      CALL DSYEV('N', 'U', N, H, N, allEig, work, lwork, tmp)
      
      if (tmp == 0) then
          print *, M,' Eigen values of H around ', E0
          call getWindow(E0, N, allEig, M, eig)
          call Reorder('A', M, eig)
          print *, eig
      else
          print *, 'Error while solving eigen value problem of H'
      end if

      end

      subroutine printEig()
      integer,parameter :: N = 100

      double precision, dimension(n)    :: allEig
      double precision, dimension(n, n) ::  H
      double precision, dimension(3*n)  :: work
      double precision, dimension(N)    :: eig
      integer :: lwork = 3*N
      integer :: tmp     
  
      common /H/  H
      

!      print *, 'Matrix'
!      print *, H

      CALL DSYEV('N', 'U', N, H, N, allEig, work, lwork, tmp)
      
      if (tmp == 0) then
          print *, ' Eigen values of H around ', allEig
      else
          print *, 'Error while solving eigen value problem of H'
      end if

      end

