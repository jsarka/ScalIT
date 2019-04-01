!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C     BLAS-LIKE subroutines for 1D vector: Complex version      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!ccccccccccccccccccccccccccccccccccccccccccccc
!c     ZEROVEC subroutine---XOUT = 0         c
!c     Initialize vector XOUT to zero        c
!ccccccccccccccccccccccccccccccccccccccccccccc
      subroutine  makeZero_CX(N, XOUT)
      implicit none
      integer,intent(IN) :: N 
      double complex,intent(OUT):: XOUT(N)  

      XOUT(1:N) = 0.0D0

      end

!ccccccccccccccccccccccccccccccccccccccccccccc
!c       Randomly generated vector           c
!ccccccccccccccccccccccccccccccccccccccccccccc
      subroutine randVec_CX(n, vec)
      implicit none
      integer, intent(IN) :: N
      double complex, intent(OUT) :: Vec(N)     

      double precision:: A(N), B(N) 

      call random_number(A)
      call random_number(B)

      Vec(1:N) = DCMPLX(A(1:N), B(1:N))     

      end 

!ccccccccccccccccccccccccccccccccccccccccc
!c     COPYVEC subroutine---XOUT = X     c
!c     Copy vector X to vector XOUT      c
!ccccccccccccccccccccccccccccccccccccccccc
      subroutine  COPY_DX(N, X, XOUT)
      implicit none
      integer,intent(IN) :: N    
      double precision, intent(IN)  :: X(N) 
      double complex, intent(OUT) :: XOUT(N) 
      
      XOUT(1:N) = X(1:N)
   
      end
!**************************************
      subroutine  COPY_CX(N, X, XOUT)
      implicit none
      integer,intent(IN) :: N    
      double complex, intent(IN)  :: X(N) 
      double complex, intent(OUT) :: XOUT(N) 
      
      XOUT(1:N) = X(1:N)

      end

!ccccccccccccccccccccccccccccccccccccccccccc
!c     AXVEC subroutine---XOUT = a X       c
!c     Rescale vector X by scalar a.       c
!ccccccccccccccccccccccccccccccccccccccccccc

      subroutine  AX_DX(N, a, X, XOUT)
      implicit none
      integer,intent(IN) :: N          
      double precision, intent(IN) :: a
      double complex, intent(IN)  :: X(N) 
      double complex, intent(OUT) :: XOUT(N) 

      XOUT(1:N) = a*X(1:N)

      end

!ccccccccccccccccccccccccccccccccccccccccccc
      subroutine  AX_CX(N, a, X, XOUT)
      implicit none
      integer,intent(IN) :: N          
      double complex, intent(IN)  :: a
      double complex, intent(IN)  :: X(N) 
      double complex, intent(OUT) :: XOUT(N) 

      XOUT(1:N) = a*X(1:N)

      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                              c 
!c     Add 2 functions to do XOUT=X+Y and XOUT=X-Y for vectors  c
!c     to reduce the time for multiplication in using AXBYVEC   c
!c     XADDY_CX(X, Y, XOUT, N)   XOUT=X+Y                       c
!c     XSUBY_CX(X, Y, XOUT, N)   XOUT=X-Y                       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
      subroutine  XADDY_CX(N, X, Y, XOUT)
      implicit none
      integer,intent(IN) :: N           
      double complex, intent(IN)  :: X(N), Y(N)
      double complex, intent(OUT) :: XOUT(N)
      
      XOUT(1:N) = X(1:N) + Y(1:N)

      end

!ccccccccccccccccccccccccccccccccccccccccccccc
!c     XSUBYVEC subroutine---XOUT = X - Y    c
!ccccccccccccccccccccccccccccccccccccccccccccc
      subroutine  XSUBY_CX(N, X, Y, XOUT)
      implicit none
      integer,intent(IN) :: N       
      double complex,intent(IN)  :: X(N), Y(N) 
      double complex,intent(OUT) :: XOUT(N)
      
      XOUT(1:N) = X(1:N) - Y(1:N)

      end

!cccccccccccccccccccccccccccccccccccccccccccccccc
!c     AXBYVEC subroutine---XOUT = a X + b Y    c
!cccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine  AXBY_DX(N, a, X, b, Y, XOUT)
      implicit none
      integer,intent(IN) :: N                  
      double precision,intent(IN) :: a, b
      double complex,intent(IN)   :: X(N), Y(N)
      double complex,intent(OUT)  :: XOUT(N) 
            
      XOUT(1:N) = a * X(1:N) + b * Y(1:N)

      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      subroutine  AXBY_CX(N, a, X, b, Y, XOUT)
      implicit none
      integer,intent(IN) :: N                  
      double complex, intent(IN) :: a, b
      double complex,intent(IN)  :: X(N), Y(N) 
      double complex,intent(OUT) :: XOUT(N) 
            
      XOUT(1:N) = a * X(1:N) + b * Y(1:N)

      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     AXBYCZVEC subroutine---XOUT = a X + b Y + c Z      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine  AXBYCZ_DX(N, a, X, b, Y, c, Z, XOUT)
      implicit none
      integer,intent(IN)           :: N   
      double precision, intent(IN) :: a, b, c            
      double complex, intent(IN)   :: X(N), Y(N), Z(N)     
      double complex, intent(OUT)  :: XOUT(N)      
             
      XOUT(1:N) = a * X(1:N) + b * Y(1:N) + c * Z(1:N)

      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine  AXBYCZ_CX(N, a, X, b, Y, c, Z, XOUT)
      implicit none
      integer,intent(IN)           :: N 
      double complex, intent(IN)   :: a, b, c              
      double complex, intent(IN)   :: X(N), Y(N), Z(N)     
      double complex, intent(OUT)  :: XOUT(N)            
      
      XOUT(1:N) = a * X(1:N) + b * Y(1:N) + c * Z(1:N)

      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     XYVEC subroutine---XOUT = X * Y (component-wise multiplication)  c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine  XY_CX(N, X, Y, XOUT)
      implicit none
      integer,intent(IN) :: N     
      double complex,intent(IN)  :: X(N), Y(N)
      double complex,intent(OUT) :: XOUT(N)

      XOUT(1:N) = X(1:N) * Y(1:N)

      end

!ccccccccccccccccccccccccccccccccccccccccccccccccc
!c     NORM function                             c
!c     Compute the norm of the vector X.         c
!ccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function NORM_CX(N, X)
      implicit none
      integer,intent(IN) ::  N                     
      double complex,intent(IN):: X(N)   

      NORM_CX = DSQRT(dble(dot_product(X(1:N),X(1:N))))

      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     NORM2 function                                    c
!c     Compute the SQUARE of the norm of the vector X.   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function NORM2_CX(N, X)
      implicit none
      integer,intent(IN)::      N      
      double complex,intent(IN):: X(N)  

      NORM2_CX = DBLE(dot_product( X(1:N), X(1:N)))

      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     DOTPROD function                                        c
!c     Compute the inner product of the two vectors X and Y.   c
!c            DOTPROD_CMPX= X^H * Y                            c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      double complex function DOTPROD_CX(N, X, Y)
      implicit none
      integer,intent(IN) :: N               
      double complex, intent(IN)  :: X(N), Y(N)
      
      DOTPROD_CX = dot_product(X(1:N) , Y(1:N))      

      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     DOTPROD function                                        c
!c     Compute the inner product of the two vectors X and Y.   c
!c         DOT_CMPX = X^T * Y                                  c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      double complex function DOT_CX(N, X, Y)
      implicit none
      integer,intent(IN) :: N      
      double complex,intent(IN)  :: X(N), Y(N)
      
      DOT_CX = sum( X(1:N) * Y(1:N))      

      end

!ccccccccccccccccccccccccccccccccccccccccccccccccc
!c             Print out vector                  c
!ccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine printVec_CX(N, vec)
      implicit none
      integer, intent(in) :: N
      double complex, intent(in) :: vec(N) 

      print *, 'Print Complex Vector: N=', N

      print *, vec(1:N)

      end subroutine

!*********************************************************
