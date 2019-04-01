!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Calculate NORM and DOTPROD in MPI environment    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function NORM_CX_MPI(COMM, NIN, X, IERR)
      implicit none
      include 'mpif.h'
      integer,intent(IN)  :: COMM, NIN
      double complex, intent(IN) :: X(NIN)
      integer,intent(OUT) :: IERR
 
      double precision :: LOC_NORM2

      LOC_NORM2 = dble(dot_product(X(1:NIN) , X(1:NIN)))

      call MPI_ALLREDUCE(LOC_NORM2,NORM_CX_MPI,1,MPI_DOUBLE_PRECISION,   &
                MPI_SUM, COMM, IERR)

      NORM_CX_MPI = DSQRT(NORM_CX_MPI)

      end

!******************************
      double precision function NORM_SVCX_MPI(COMM, NIN, X, IERR)
      implicit none
      integer,intent(IN)  :: COMM, NIN
      double complex, intent(IN) :: X(NIN)
      integer,intent(OUT) :: IERR      
 
      call NORM2_MVCX_MPI(COMM, NIN,1, X, NORM_SVCX_MPI, IERR)

      NORM_SVCX_MPI = DSQRT(NORM_SVCX_MPI)    

      end

!*******************************************
      subroutine NORM_MVCX_MPI(COMM, NIN, NUM, X, NORM_MPI, IERR)
      implicit none
      integer,intent(IN)  :: COMM, NIN, NUM
      double complex,intent(IN) :: X(NIN,NUM)
      double precision,intent(OUT) :: NORM_MPI(NUM)
      integer,intent(OUT) :: IERR      
 
      call NORM2_MVCX_MPI(COMM, NIN, NUM, X, NORM_MPI,IERR)

      NORM_MPI(1:NUM) = DSQRT(NORM_MPI(1:NUM))    

      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                   c
!c  Calculate NORM and DOTPROD in MPI environment    c
!c                                                   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function NORM2_CX_MPI(COMM, NIN, X, IERR)
      implicit none
      include 'mpif.h'
      integer,intent(IN)  :: COMM, NIN
      double complex, intent(IN) :: X(NIN)
      integer,intent(OUT) :: IERR

      double precision :: LOC_NORM2

      LOC_NORM2 = dble(dot_product(X(1:NIN) , X(1:NIN)))

      call MPI_ALLREDUCE(LOC_NORM2,NORM2_CX_MPI,1,MPI_DOUBLE_PRECISION,   &
                MPI_SUM, COMM, IERR)
      end


!****************************
      double precision function NORM2_SVCX_MPI(comm, NIN, X, IERR)
      implicit none    
      integer,intent(IN)  :: COMM, NIN
      double complex, intent(IN) :: X(NIN)
      integer,intent(OUT) :: IERR

      call NORM2_MVCX_MPI(COMM, NIN, 1, X, NORM2_SVCX_MPI, IERR)

      end
!****************************
      subroutine NORM2_MVCX_MPI(COMM, NIN, NUM, X, NORM2, IERR)
      implicit none
      include 'mpif.h'
      integer,intent(IN)  :: COMM, NIN, NUM
      double complex, intent(IN) :: X(NIN,NUM)
      double precision, intent(OUT) :: NORM2(NUM)
      integer,intent(OUT) :: IERR

      double precision :: LOC_NORM2(NUM)
      integer :: I

      do I = 1, NUM
           LOC_NORM2(I) = dble(dot_product(X(1:NIN, I) , X(1:NIN, I)))
      end do

      call MPI_ALLREDUCE(LOC_NORM2, NORM2, NUM,MPI_DOUBLE_PRECISION,   &
                MPI_SUM, COMM, IERR)      
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Calculate NORM and DOTPROD in MPI environment     c
!c         DOTPROD = X^H * Y                          c  
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double complex function DOTPROD_CX_MPI(COMM, NIN, X, Y, IERR)
      implicit none
      include 'mpif.h'
      integer,intent(IN)  :: COMM, NIN
      double complex, intent(IN) :: X(NIN), Y(NIN)
      integer,intent(OUT) :: IERR

      double complex :: LOC_DOT

      LOC_DOT = dot_product(X(1:NIN) , Y(1:NIN))

      call MPI_ALLREDUCE(LOC_DOT,DOTPROD_CX_MPI,1,MPI_DOUBLE_COMPLEX,   &
                MPI_SUM, COMM, IERR)

      end

!ccccccccccccccccccccccccc
      double complex function DOTPROD_SVCX_MPI(COMM, NIN, X, Y, IERR)
      implicit none    
      integer,intent(IN)  :: COMM, NIN
      double complex, intent(IN) :: X(NIN), Y(NIN)
      integer,intent(OUT) :: IERR

      call DOTPROD_MVCX_MPI(COMM,NIN,1,X,1,Y,DOTPROD_SVCX_MPI,IERR)

      end

!***********************
      subroutine DOTPROD_MVCX_MPI(COMM,NIN,NUMX,X,NUMY,Y,DOTCX,IERR)
      implicit none
      include 'mpif.h'
      integer,intent(IN)  :: COMM, NIN, NUMX, NUMY
      double complex, intent(IN) :: X(NIN, NUMX),Y(NIN, NUMY)
      double complex, intent(OUT):: DOTCX(NUMX, NUMY) 
      integer,intent(OUT) :: IERR

      double complex :: LOC_DOT(numx, numy) 
      integer :: I, J

      do I = 1, NUMX
          do J = 1, NUMY
             LOC_DOT(I, J) = dot_product(X(1:NIN, I) , Y(1:NIN, J))
         end do
      end do

      call MPI_ALLREDUCE(LOC_DOT,DOTCX, numx*numy,MPI_DOUBLE_COMPLEX,   &
                MPI_SUM, COMM, IERR)

      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Calculate NORM and DOTPROD in MPI environment     c
!c         DOTPROD = X^T * Y                          c  
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double complex function DOT_CX_MPI(COMM, NIN, X, Y, IERR)
      implicit none
      include 'mpif.h'
      integer,intent(IN)  :: COMM, NIN
      double complex, intent(IN) :: X(NIN), Y(NIN)
      integer,intent(OUT) :: IERR

      double complex :: LOC_DOT

      LOC_DOT = sum( X(1:NIN) * Y(1:NIN))

      call MPI_ALLREDUCE(LOC_DOT,DOT_CX_MPI,1,MPI_DOUBLE_COMPLEX,   &
                MPI_SUM, COMM, IERR)

      end

!******************************
      double complex function DOT_SVCX_MPI(COMM, NIN, X, Y, IERR)
      implicit none    
      integer,intent(IN)  :: COMM, NIN
      double complex, intent(IN) :: X(NIN), Y(NIN)
      integer,intent(OUT) :: IERR
 
      call DOT_MVCX_MPI(COMM, NIN, 1, X, 1, Y, DOT_SVCX_MPI, IERR)

      end

!****************************
      subroutine DOT_MVCX_MPI(COMM, NIN, NUMX, X, NUMY, Y, DOTCX, IERR)
      implicit none
      include 'mpif.h'
      integer,intent(IN)  :: COMM, NIN, NUMX, NUMY
      double complex, intent(IN) :: X(NIN, NUMX),Y(NIN, NUMY)
      double complex, intent(OUT):: DOTCX(NUMX, NUMY) 
      integer,intent(OUT) :: IERR

      double complex :: LOC_DOT(numx, numy) 
      integer        :: i, j

      do I =1 , NUMX
           do J = 1, NUMY
              LOC_DOT = sum( X(1:NIN, I) * Y(1:NIN, J))
           end do
      end do

      call MPI_ALLREDUCE(LOC_DOT,DOTCX, NUMX*numY,MPI_DOUBLE_COMPLEX,   &
                MPI_SUM, COMM, IERR)

      end

!******************************************
