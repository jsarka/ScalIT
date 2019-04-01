!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C     BLAS-LIKE subroutines for 1D vector: Real Version         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!ccccccccccccccccccccccccccccccccccccccccccccc
!c     ZEROVEC subroutine---XOUT = 0         c
!c     Initialize vector XOUT to zero        c
!ccccccccccccccccccccccccccccccccccccccccccccc
   subroutine  makeZero(N, XOUT)
      implicit none
      integer,intent(IN) :: N                          
      double precision,intent(OUT):: XOUT(N)

      XOUT(1:N) = 0.0d0

      end

!ccccccccccccccccccccccccccccccccccccccccccc
!c       Randomly generated vector         c
!ccccccccccccccccccccccccccccccccccccccccccc
      subroutine randVec(n, vec)
      implicit none
      integer, intent(IN) :: N
      double precision, intent(OUT) :: Vec(N)     

      call random_number(vec)

      end 

!ccccccccccccccccccccccccccccccccccccccccc
!c     COPYVEC subroutine---XOUT = X     c
!c     Copy vector X to vector XOUT      c
!ccccccccccccccccccccccccccccccccccccccccc
      subroutine  COPY(N, X, XOUT)
      implicit none
      integer,intent(IN) :: N                
      double precision,intent(IN)  :: X(N)
      double precision,intent(OUT) :: XOUT(N)
      
      XOUT(1:N) = X(1:N)

      end

!ccccccccccccccccccccccccccccccccccccccccccc
!c     AXVEC subroutine---XOUT = a X       c
!c     Rescale vector X by scalar a.       c
!ccccccccccccccccccccccccccccccccccccccccccc
      subroutine  AX(N, A, X, XOUT)
      implicit none
      integer,intent(IN) :: N     
      double precision, intent(IN) :: A    
      double precision, intent(IN) :: X(N)
      double precision, intent(OUT):: XOUT(N)

      XOUT(1:N) = A*X(1:N)

      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Add 2 functions to do XOUT=X+Y and XOUT=X-Y for vectors  c
!c     to reduce the time for multiplication in using AXBYVEC   c
!c     XADDYVEC(X, Y, XOUT, N)   XOUT=X+Y                       c
!c     XSUBYVEC(X, Y, XOUT, N)   XOUT=X-Y                       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
      subroutine  XADDY(N, X, Y, XOUT)
      implicit none
      integer,intent(IN) :: N    
      double precision, intent(IN)  :: X(N), Y(N)
      double precision, intent(OUT) :: XOUT(N)
      
      XOUT(1:N) = X(1:N) + Y(1:N)

      end

!ccccccccccccccccccccccccccccccccccccccccccccc
!c     XSUBYVEC subroutine---XOUT = X - Y    c
!ccccccccccccccccccccccccccccccccccccccccccccc
      subroutine  XSUBY(N, X, Y, XOUT)
      implicit none
      integer,intent(IN) :: N                
      double precision, intent(IN)    :: X(N), Y(N)
      double precision, intent(OUT)   :: XOUT(N)
      
      XOUT(1:N) = X(1:N) - Y(1:N)

      end

!cccccccccccccccccccccccccccccccccccccccccccccccc
!c     AXBYVEC subroutine---XOUT = a X + b Y    c
!cccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine  AXBY(N, A, X, B, Y, XOUT)
      implicit none
      integer,intent(IN) :: N               
      double precision,intent(IN)   :: A, B    
      double precision, intent(IN)  :: X(N), Y(N)
      double precision, intent(OUT) :: XOUT(N)
            
      XOUT(1:N) = A * X(1:N) + B * Y(1:N)

      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     AXBYCZVEC subroutine---XOUT = a X + b Y + c Z      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine  AXBYCZ(N, A, X, B, Y, C, Z, XOUT)
      implicit none
      integer,intent(IN) :: N    
      double precision, intent(IN)   :: a, b, c             
      double precision, intent(IN)   :: X(N), Y(N), Z(N)                                 
      double precision, intent(OUT)  :: XOUT(N) 
      
      XOUT(1:N) = a * X(1:N) + b * Y(1:N) + c * Z(1:N)

      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     XYVEC subroutine---XOUT = X * Y (component-wise multiplication)  c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine  XY(N, X, Y, XOUT)
      implicit none
      integer,intent(IN) :: N  
      double precision, intent(IN)  :: X(N), Y(N)
      double precision, intent(OUT) :: XOUT(N)

      XOUT(1:N) = X(1:N) * Y(1:N)

      end

!ccccccccccccccccccccccccccccccccccccccccccccccccc
!c     NORM function                             c
!c     Compute the norm of the vector X.         c
!ccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function NORM(N, X)
      implicit none
      integer,intent(IN) ::  N                      
      double precision,intent(IN):: X(N)  

      NORM = DSQRT(dot_product(X(1:N), X(1:N)))

      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     NORM2 function                                    c
!c     Compute the SQUARE of the norm of the vector X.   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function NORM2(N, X)
      implicit none
      integer,intent(IN)::      N             
      double precision,intent(IN):: X(N)   

      NORM2 = dot_product( X(1:N), X(1:N))

      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     DOTPROD function                                        c
!c     Compute the inner product of the two vectors X and Y.   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function DOTPROD(N, X, Y)
      implicit none
      integer,intent(IN) :: N                   
      double precision, intent(IN)  :: X(N), Y(N)
      
      DOTPROD = dot_product(X(1:N) , Y(1:N))      

      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     DOTPROD function                                        c
!c     Compute the result of Vi^T * D * Vj                     c
!c    Vi, Vj are vectors, and D is the diagonal matrix         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      double precision function DOTVTDV(N, Vi, D, Vj)
      implicit none
      integer,intent(IN) :: N               
      double precision,intent(IN)  :: Vi(N), D(N), Vj(N)
      
      DOTVTDV = sum(Vi(1:N)*D(1:N)*Vj(1:N))

      end

!cccccccccccccccccccccccccccccccccccccccccccc
!c         Print the vector                 c
!cccccccccccccccccccccccccccccccccccccccccccc

      subroutine printVec(N, vec)
      implicit none
      integer, intent(in) :: N
      double precision, intent(in) :: vec(N)

      print *, 'Print Real Vector: N=', N

      print *, vec(1:N)

      end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
