!ccccccccccccccccccccccccccccccccccccccc
!c  
!c        testing program for qmr
!c
!ccccccccccccccccccccccccccccccccccccccc

      program qmr_drc

      integer,parameter :: N = 23
      integer, parameter :: QMRMAX = 1000
      double precision, parameter :: QMRTOL = 1.0D-4

      double complex, dimension(N):: b, x0, x   
      
      integer   i,iter, qmr  
      double precision res
      external :: HMULT, PMULT

      print *, ' test for qmr subroutine'
      

      call HINIT()

      do i=1,n
         x0(i) = i * 1.0d0 
      end do

      print *, 'Calculate and Display X, HX, PX'

      call HMULT(23, X0,B)
      call PMULT(23, X0,X)

      do i = 1, n
         print *, x0(i), b(i), X(i)
      end do
      
      print *
      print *, ' Now performing QMR'
      iter =  QMR_CX(0, QMRMAX,QMRTOL,N,B,HMULT, PMULT, X,res)
      print *, 'Iteration number: ',iter
      print *, 'RES: ', res

      print *
      print *, 'Initial vector x0, and solution x, relative error'
      do i=1,n
         print *, x0(i), x(i), (x0(i)-x(i))/x0(i)
         print *
      end do 
      
      print *, 'Finish QMR testing'

      END


!ccccccccccccccccccccccccccccccccccccccc
!c  
!c        testing functions for qmr
!c
!ccccccccccccccccccccccccccccccccccccccc

!c
!c     Get symmetric h matrix
!c
      SUBROUTINE  HINIT()
      
      integer,parameter :: N  = 23 
      double complex, dimension(N, N):: H                    

      common /H/  H

      integer i, j
      
      do I=1,N
         do J = 1, I-1
            H(I, J)  = rand()
            H(J, I)  = H(I, J)
         end do
	H(I, I) = DCMPLX(rand(), rand()) 
      end do    

      end

!c
!c     y = H x
!c
      SUBROUTINE  HMULT(M, X,Y)
      
      integer,parameter :: N = 23

      integer :: M
      double complex, dimension(n)    :: X, Y
      double complex, dimension(n, n) ::  H

      common /H/ H
    
      integer i, j

!c      print *, 'call HMult'

      do i = 1, N
         y(i) = 0.0d0
         do j=1,N
            y(i) = y(i) + H(i,j)*x(j)
         end do
      end do 
         
      end

!c
!c     y = P^(-1) x, Jacobi preconditioning
!c
      SUBROUTINE  PMULT(M, X,Y)
   
      integer,parameter :: N = 23

      integer :: M

      double complex, dimension(n)    :: X, Y
      double complex, dimension(n, n) ::  H

      common /H/  H
    
      integer i
      
!c      print *, 'Call PMult'

      do  i = 1, N
         y(i) = x(i)/H(i,i)
      end do
         
      end

