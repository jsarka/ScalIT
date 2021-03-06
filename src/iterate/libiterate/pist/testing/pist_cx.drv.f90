!ccccccccccccccccccccccccccccccccccccccc
!c  
!c        testing program for PIST_CX
!c
!ccccccccccccccccccccccccccccccccccccccc

      program pist_cx_drv
      implicit none       
      integer,parameter  :: N = 100
      integer, parameter :: EIGNUM = 4

      double complex, dimension(EIGNUM) ::  EIG
      double complex, dimension(N)      ::  x   
      external :: hmult, pmult, qmrmult    

      double precision, parameter :: E0     = 0.0D0
      double precision, parameter :: ETOL   = 1.0e-4


      integer  ::  i,num
      double precision :: res
      integer  :: testType = 0;
      integer  :: MAX_NUM  = 50
      integer  :: pist1_cx, pist2_cx, pist51_cx, pist52_cx

      print *, ' test for Complex version of PIST subroutine'
      print *, ' Input testing type:'
      print *, ' -1/0: Use ZGEEV directly'
      print *, ' 1/2: Faster PISTCONV'
      print *, ' 3/4: Slower PISTCONV'
      read(*, *) testType
 
      call HINIT()
!      res = max(X(1:N))

      do i=1,n
         x(i) = i * 1.0d0 
      end do

     select case (testType)
     case (-1)
          print *
          print *, '  Get Part of Eigen values directly'
          print *, ' Input number of eigen values:'      
          read(*,*) num
          call printEig0_CX(E0, num)


     case (0)
      print *
      print *, ' Get the eigen value directly using ZGEEV'
      call printEig_CX()

     case (1)
      print *
      print *, '  Performing PIST using Faster Convergent Criterion'
      if(pist1_cx(E0, ETOL, N, X, EIGNUM, MAX_NUM, HMULT,QMRMULT, EIG, RES)>0) THEN
          print *, 'Res:', res
          print *, 'Eig:', EIG
      else
         print *, 'Error in PISTCONV'
      end if

     case(2)
      print *
      print *, ' Performing PIST using Faster Convergent Criterion'
      if(pist2_cx(E0, ETOL, N, X, EIGNUM, MAX_NUM, HMULT, QMRMULT,EIG, RES)>0) THEN
         print *, 'Res:', res
         print *, 'Eig: ', eig
      else
         print *, 'Error in PISTCONV'
      end if

    case (3)
      print *
      print *, ' Performing PIST using Slower Convergent Criterion'
      if(pist51_cx(E0, ETOL, N, X, EIGNUM, MAX_NUM, HMULT, QMRMULT,EIG, RES)>0) THEN
         print *, 'Res:', res
         print *, 'Eig: ', eig
      else
         print *, 'Error in PISTCONV'
      end if

    case (4:)
      print *
      print *, ' Performing PIST using Slower Convergent Criterion'
      if(pist52_cx(E0, ETOL, N, X, EIGNUM, MAX_NUM, HMULT, QMRMULT, EIG, RES)>0) THEN
         print *, 'Res:', res
         print *, 'Eig: ', eig
      else
         print *, 'Error in PISTCONV'
      end if
    end select 

      print *, 'Finish PIST testing'

      END

!c
!c     Get symmetric h matrix
!c
      SUBROUTINE  HINIT()
      
      integer,parameter :: N  = 100 
      double complex, dimension(N, N):: H                    

      common /H/  H

      integer i, j
      
      do I=1,N
         do J = 1, I
            H(I, J)  = rand()
            H(J, I)  = H(I, J)
         end do
      end do    

      do I = 1, N
          H(I, I) = DCMPLX(rand(), rand())
      end do

      end

!c
!c     y = H x
!c
      SUBROUTINE  HMULT(m, X,Y)
      
      integer,parameter :: N = 100

      double complex, dimension(n)    :: X, Y
      double complex, dimension(n, n) ::  H

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
      SUBROUTINE  PMULT(m, X,Y)

      integer,parameter :: N = 100

      double complex, dimension(n)    :: X, Y
      double complex, dimension(n, n) ::  H

      common /H/  H
    
      integer i
      
!c      print *, 'Call PMult'

      do  i = 1, N
         y(i) = x(i)/DBLE(H(i,i))
      end do
         
      end

      integer function QMRMULT(M, B, X, RES)
      implicit none
      integer, parameter :: QMRMAX   = 1000
      integer, parameter :: QMR_TYPE = 0
      integer, parameter :: N = 100
      double precision, parameter :: QMRTOL = 1.0D-4      

      integer, intent(IN)  :: M
      double complex, dimension(M), intent(IN)  :: B
      double complex, dimension(M), intent(OUT) :: X
      double precision, intent(OUT) :: RES

      external :: HMULT, PMULT
      integer, external :: QMR_CX
      
      QMRMULT = QMR_CX(QMR_TYPE, QMRMAX,QMRTOL, N, B, HMULT, PMULT, X, RES)

      end

!c
!c   Get the eigen values of H directly using ZSYEV
!c
!      subroutine printEig_CX1()
!      integer,parameter :: N = 100
!      integer, parameter :: LDVL=1, LDVR=1
!      double complex :: VL, VR
!      double complex, dimension(n)    :: eig
!      double complex, dimension(n, n) ::  H
!      double complex, dimension(3*n)  :: work
!      double complex, dimension(3*n)  :: RWORK
!      integer :: lwork = 3*N
!      integer :: tmp
  
!      common /H/  H
      
      ! CALL DSYEV('N', 'U', N, H, N, Eig, work, lwork, tmp)
!      CALL ZSYEV('N', 'U', N, H, N, Eig, WORK, LWORK, RWORK, tmp)
      
!      if (tmp == 0) then
!          print *, 'Eigen values of H using ZDYSV'
!          print *, eig
!      else
!          print *, 'Error while solving eigen value problem of H'
!      end if

!      end


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
          call Reorder_CX('A', 1, Num, eig0)
          print *, eig0
      else
          print *, 'Error while solving eigen value problem of H'
      end if

      end
