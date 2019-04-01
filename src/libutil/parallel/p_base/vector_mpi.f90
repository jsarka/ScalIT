!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Calculate NORM and DOTPROD in MPI environment    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccc
!c    norm for One Vector      c
!ccccccccccccccccccccccccccccccc
      double precision function NORM_MPI(COMM, NIN, X, IERR)
      implicit none
      include 'mpif.h' 
      integer,intent(IN)  :: COMM, NIN
      double precision, intent(IN) :: X(NIN)
      integer,intent(OUT) :: IERR

      double precision :: LOC_NORM2

      LOC_NORM2 = dot_product(X(1:NIN) , X(1:NIN))

      call MPI_ALLREDUCE(LOC_NORM2,NORM_MPI,1,MPI_DOUBLE_PRECISION,   &
                MPI_SUM, COMM, IERR)

      NORM_MPI = DSQRT(NORM_MPI)

      end

!ccccccccccccccccccccccccccccccc
!c    norm for Many Vectors    c
!ccccccccccccccccccccccccccccccc
      double precision function NORM_SV_MPI(COMM, NIN, X, IERR)
      implicit none
      integer,intent(IN)  :: COMM, NIN
      double precision, intent(IN) :: X(NIN)
      integer,intent(OUT) :: IERR     
 
      call NORM2_MV_MPI(COMM, NIN, 1, X, NORM_SV_MPI, IERR)

      NORM_SV_MPI = DSQRT(NORM_SV_MPI)

      end
!*********************
      subroutine NORM_MV_MPI(COMM, NIN, NUM, X, NORM, IERR)
      implicit none
      integer,intent(IN)  :: COMM, NIN, NUM
      double precision, intent(IN) :: X(NIN,NUM)
      double precision, intent(out):: norm(NUM)
      integer,intent(OUT) :: IERR     
 
      call NORM2_MV_MPI(COMM, NIN, NUM, X, NORM, IERR)

      NORM(1:NUM) = DSQRT(NORM(1:NUM))

      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                   c
!c  Calculate NORM and DOTPROD in MPI environment    c
!c                                                   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function NORM2_MPI(COMM, NIN, X, IERR)
      implicit none
      include 'mpif.h'
      integer,intent(IN)  :: COMM, NIN
      double precision, intent(IN) :: X(NIN)
      integer,intent(OUT) :: IERR

      double precision :: LOC_NORM2

      LOC_NORM2 = dot_product(X(1:NIN) , X(1:NIN))

      call MPI_ALLREDUCE(LOC_NORM2,NORM2_MPI,1,MPI_DOUBLE_PRECISION,   &
                MPI_SUM, COMM, IERR)

      end

!************************************
      double precision function NORM2_SV_MPI(COMM, NIN, X, IERR)
      implicit none      
      integer,intent(IN)  :: COMM, NIN
      double precision, intent(IN) :: X(NIN)
      integer,intent(OUT) :: IERR
 
      call NORM2_MV_MPI(COMM, NIN, 1, X, NORM2_SV_MPI, IERR) 
      
      end

!************************************
      subroutine NORM2_MV_MPI(COMM, NIN, NUM, X, NORM2, IERR)
      implicit none
      include 'mpif.h'
      integer,intent(IN)  :: COMM, NIN, NUM
      double precision, intent(IN)  :: X(NIN,NUM) 
      double precision, intent(OUT) :: NORM2(NUM) 
      integer,intent(OUT) :: IERR

      double precision :: LOC_NORM2(NUM)
      integer :: I

      do I = 1, NUM
          LOC_NORM2(I) = dot_product(X(1:NIN, I) , X(1:NIN, I))
      end do

      call MPI_ALLREDUCE(LOC_NORM2,NORM2, NUM, MPI_DOUBLE_PRECISION,   &
                MPI_SUM, COMM, IERR)

      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Calculate NORM and DOTPROD in MPI environment     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function DOTPROD_MPI(COMM, NIN, X, Y, IERR)
      implicit none
      include 'mpif.h'
      integer,intent(IN)  :: COMM, NIN
      double precision, intent(IN) :: X(NIN), Y(NIN)
      integer,intent(OUT) :: IERR

      double precision :: LOC_DOT

      LOC_DOT = dot_product(X(1:NIN) , Y(1:NIN))

      call MPI_ALLREDUCE(LOC_DOT,DOTPROD_MPI,1,MPI_DOUBLE_PRECISION,   &
                MPI_SUM, COMM, IERR)
      end

!****************************
      double precision function DOTPROD_SV_MPI(COMM, NIN, X, Y, IERR)
      implicit none      
      integer,intent(IN)  :: COMM, NIN
      double precision, intent(IN) :: X(NIN), Y(NIN)  
      integer,intent(OUT) :: IERR

      call DOTPROD_MV_MPI(COMM, NIN, 1, X, 1, Y, DOTPROD_SV_MPI, IERR)

      end

!*****************************
      subroutine DOTPROD_MV_MPI(COMM,NIN,NUMX,X,NUMY,Y,DOTMAT,IERR)
      implicit none
      include 'mpif.h'
      integer,intent(IN)  :: COMM, NIN, NUMX, NUMY
      double precision, intent(IN)  :: X(NIN,NUMX), Y(NIN,NUMY)
      double precision, intent(OUT) :: DOTMAT(NUMX, NUMY)
      integer,intent(OUT) :: IERR

      double precision :: LOC_DOT(NUMX, NUMY)
      integer :: I, J

      do I = 1, NUMX
           do J = 1, NUMY
                LOC_DOT(I, J) = dot_product(X(1:NIN, I) , Y(1:NIN, J))
           end do
      end do

      call MPI_ALLREDUCE(LOC_DOT,DOTMAT, NUMX*NUMY,MPI_DOUBLE_PRECISION,   &
                MPI_SUM, COMM, IERR)
      end
!*****************************************************************************
