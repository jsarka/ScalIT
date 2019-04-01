!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                   c
!c   testing program for LAN for Hermitian matrix    c
!c                                                   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc 

      program lanczos_dx_drv

      integer,parameter  :: N = 100
      integer, parameter :: EIGNUM = 4
      double precision, parameter :: E0     = -10.0D0
      double precision, parameter :: ETOL   = 1.0e-7
      double precision :: ct1, ct2

      double precision, dimension(EIGNUM) ::  EIG
      double complex, dimension(N)      ::  x   
      external :: hmult_cx
      integer  :: lanczos0, lanczos1, lanczos2, lanczos3

      integer  ::  i,iter, num  
      double precision :: res
      integer  :: testType = 0;
      integer  :: MAX_NUM  = 100

      print *, ' test for LANCZOS subroutine for Hermitian Matrix'
      print *, ' Input testing type: '
      print *, ' -1/0: for DSYEV directly'
      print *, ' 1,2:  for fast LANCZOSCONV '
      print *, ' 3,4:  for slow LANCZOSCONV '     
      read(*, *) testType
 
      call HINIT_CX()
!      res = max(X(1:N))

      do i=1,n
         x(i) = i * 1.0d0 
      end do

     call cpu_time(ct1)
     select case (testType)
    case (-1)
      print *
      print *, '  Get Part of Eigen values directly'
!      print *, ' Input number of eigen values:'      
!      read *, num
      call printEig0_CX(E0, eignum)

    case (0)
      print *
      print *, 'Get Eigen values directly:'      
      call printEig_cx()

     case(1)
     print *
      print *, ' Performing LANCZOS using faster Convergent Criterion'
      if(lanconv1_dx(E0, ETOL, N, X, EIGNUM, MAX_NUM, HMULT_cx, EIG, RES)>0) THEN
         print *, 'Res:', res
         print *, 'Eig: ', eig
      else
         print *, 'Error in LANCONV1'
         print *, 'Res:', res
         print *, 'Eig: ', eig
      end if

     case(2)
     print *
      print *, ' Performing LANCZOS using Faster Convergent Criterion'
      if(lanconv2_dx(E0, ETOL, N, X, EIGNUM, MAX_NUM, HMULT_cx, EIG, RES)>0) THEN
         print *, 'Res:', res
         print *, 'Eig: ', eig
      else
         print *, 'Error in LANCONV2'
         print *, 'Res:', res
         print *, 'Eig: ', eig
      end if

     case(3)
     print *
      print *, ' Performing less strict LANCZOS using Slower Convergent Criterion'
      if(lanConv51_dx(E0, ETOL, N, X, EIGNUM, MAX_NUM,HMULT_cx, EIG, RES)>0) THEN
         print *, 'Res:', res
         print *, 'Eig: ', eig
      else
         print *, 'Error in LANCZOSCONVERG'
         print *, 'Res:', res
         print *, 'Eig: ', eig
      end if


    case (4)
      print *
      print *, ' Performing less strict LANCZOS using Slower Convergent Criterion'
      if(lanconv52_dx(E0, ETOL, N, X, EIGNUM, MAX_NUM, HMULT_cx, EIG, RES)>0) THEN
         print *, 'Res:', res
         print *, 'Eig: ', eig
      else
         print *, 'Error in LANCZOSCONVERG'
         print *, 'Res:', res
         print *, 'Eig: ', eig
      end if
    end select 

    call cpu_TIME(ct2)

    print *, 'CPU Time:', ct2-ct1
    print *, 'Finish LANCZOS testing'

    END

!c
!c     Get Hermitian H matrix
!c
      SUBROUTINE  HINIT_cx()
      
      integer,parameter :: N  = 100 
      double complex, dimension(N, N):: H                    

      common /H/  H

      integer i, j
      
      do I=1,N
         do J = 1, I
            H(I, J)  = DCMPLX(rand(),rand())
            H(J, I)  = DCONJG(H(I, J))
         end do
      end do    

      do I=1,N
          H(I,I) = rand()
      end do

      end

!c
!c     y = H x
!c
      SUBROUTINE  HMULT_CX(m, X,Y)
      
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
      SUBROUTINE  PMULT_cx(m, X,Y)

      integer,parameter :: N = 100

      double complex, dimension(n)    :: X, Y
      double complex, dimension(n, n) ::  H

      common /H/  H
    
      integer i
      
!c      print *, 'Call PMult'

      do  i = 1, N
         y(i) = x(i)/H(i,i)
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

