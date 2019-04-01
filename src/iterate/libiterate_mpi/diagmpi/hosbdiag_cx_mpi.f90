!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C     HOSBDIAG_MPI ROUTINES                                                   C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 subroutine HOSBDIAG_CX_MPI(COMM, BJMAX, BJTOL, BLK, MAK, NDK,   &
                         HOSB, VOSB, RES, IERR)
      implicit none
      include 'mpif.h'
      integer, intent(IN) ::  COMM, BJMAX, BLK, MAK, NDK
      double precision, intent(IN) :: BJTOL
      double complex,intent(INOUT) :: HOSB(NDK, MAK, MAK, BLK)
      double precision,intent(OUT) :: VOSB(MAK*MAK, BLK)
      double complex,intent(OUT)   :: RES(NDK, MAK, BLK)
      integer,  intent(OUT) ::  IERR    

!ccccccccccccccccccccccccccccccccc
      double precision  :: YE, LOCYE, INITYE   !error as the converg. criterior
      double precision  :: YERROR, DYE         !error as the converg. criterior

      double precision  ::  sumOffDiag_CX_MPI
      double precision  ::  PHI, MJACOBI_CX_MPI

      integer :: ITER, I, J, K   

      MAIN_LOOP : do ITER = 1, BLK 

         INITYE = sumOffDiag_CX_MPI(COMM, NDK,MAK,HOSB(1,1,1,ITER), IERR)         
         
         YE = INITYE;   DYE = INITYE;      K = 1

         call VINIT(MAK, VOSB(1, ITER))
     
         do 10 while((DABS(DYE).gt.INITYE*BJTOL).and.(K.le.BJMAX))

            do I = 1,(MAK-1)   
               do J = 1,(MAK-I)     
                  PHI=MJACOBI_CX_MPI(COMM,J,J+I,MAK,NDK, HOSB(1,1,1,iter),IERR)
                  call VUPDATE(J, J+I, PHI, MAK, VOSB(1,ITER))                  
               end do
            end do              ! One Jacobi cycle                      

            DYE = sumOffDiag_CX_MPI(COMM,NDK,MAK,HOSB(1,1,1,ITER), IERR) - YE           
            YE = YE + DYE;           K = K+1

 10      end do            

         do I=1, MAK  
            RES(1:NDK,I,ITER) = HOSB(1:NDK,I,I,ITER)
         end do
    
         do I = 1, MAK
            HOSB(1:NDK, i, i, ITER) = 0.0D0
         end do

      end do MAIN_LOOP 
                       
 end
!****************************************************************************

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                      YERROR function:                                    c
!c     Calculates the total error in the off-block-diagonal elements of H   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  double precision function sumOffDiag_CX_MPI(COMM, NIN,NOUT,HOUT, IERR)
      implicit none      
      include 'mpif.h'
      integer, intent(IN)    :: COMM,NIN, NOUT             
      double complex,intent(IN) :: HOUT(NIN,NOUT,NOUT)
      integer, intent(OUT)   :: IERR

      integer  :: I, J       
      double precision :: locSum

      locSum = 0.0D0
      do I = 1, NOUT
         do J = 1,I-1
            locSum = locSum + abs(dot_product(HOUT(1:NIN,J,I),HOUT(1:NIN,J,I)))
         end do
      end do

      call MPI_ALLREDUCE(locSum, sumOffDiag_CX_MPI,1,MPI_DOUBLE_PRECISION, &
             MPI_SUM, COMM, ierr)

      sumOffDiag_CX_MPI = 2.0d0 * sumOffDiag_CX_MPI

  end
!********************************************************************

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  double precision function sumDiag_CX_MPI(COMM, NIN,NOUT,HOUT, IERR)
      implicit none      
      include 'mpif.h'
      integer, intent(IN)    :: COMM,NIN, NOUT             
      double complex,intent(IN) :: HOUT(NIN,NOUT,NOUT)
      integer, intent(OUT)   :: IERR

      integer  :: I
      double precision :: locSum

      locSum = 0.0D0
      do I = 1, NOUT   
          locSum = locSum + abs(dot_product(HOUT(1:NIN,I,I),HOUT(1:NIN,I,I)))        
      end do

      call MPI_ALLREDUCE(locSum, sumDiag_CX_MPI,1,MPI_DOUBLE_PRECISION, &
             MPI_SUM, COMM, ierr)
  end
!********************************************************************

