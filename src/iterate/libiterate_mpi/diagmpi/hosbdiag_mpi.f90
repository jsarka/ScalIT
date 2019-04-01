!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C     HOSBDIAG_MPI ROUTINES                                                   C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine HOSBDIAG_MPI(COMM, BJMAX, BJTOL, BLK, MAK, NDK,   &
                          HOSB, VOSB, RES, IERR, ID)
      implicit none
      include 'mpif.h'

      integer, intent(IN) :: COMM, BJMAX, Blk, Mak, NDK, ID
      double precision, intent(IN) :: BJTOL
      double precision, intent(INOUT) :: HOSB(NDK, MAK, MAK, BLK)
      double precision, intent(OUT) :: VOSB(MAK*MAK, BLK),RES(NDK, MAK, BLK)
      integer,  intent(OUT) :: IERR 

!ccccccccccccccccccccccccccccccccc
      double precision :: YE, LOCYE, INITYE   !error as the converg. criterior
      double precision :: YERROR, DYE         !error as the converg. criterior
      double precision :: eff(5), eff0

      double precision :: sumOffDiag_MPI, PHI, MJACOBI_MPI, sumDiag_MPI 
      integer :: ITER, I, J, K, rootID  
      double precision :: ct1s, ct2s, ct1a, ct2a, ct3a, ct4a!, MPI_WTime

      if (id==rootID) print *

      rootID = 0
!ccccccccccccccc  Start Coding ccccccccccccccccccccccccc
      MAIN_LOOP : do ITER = 1, BLK        ! NUMBER OF JACOBI TRANSFORMS   

         INITYE = sumOffDiag_MPI(COMM, NDK,MAK,HOSB(1,1,1,ITER), IERR)         
         
         YE = INITYE;   DYE = INITYE;  K = 1

         call VINIT(MAK, VOSB(1, ITER))
     
         if (id==rootID) ct1s = MPI_WTIME()
         !if (id==rootID) print *

         do 10 while((DABS(DYE) > INITYE*BJTOL).and.(K <= BJMAX))
         if (id==rootID) ct1a = MPI_WTIME()

            do I = 1,(MAK-1)   
               do J = 1,(MAK-I)     
                  PHI = MJACOBI_MPI(COMM,J,J+I,MAK,NDK, HOSB(1,1,1,iter),IERR)
                  call VUPDATE(J, J+I, PHI, MAK, VOSB(1,ITER))                  
               end do
            end do              ! One Jacobi cycle                      
            if (id==rootID) ct2a = MPI_WTIME()
            !if (id==rootID) print *, 'MPI WTime MPI:', K, ct2a-ct1a

            DYE = sumOffDiag_MPI(COMM,NDK,MAK,HOSB(1,1,1,ITER),IERR) - YE           
            YE = YE + DYE;    K = K+1

            if ((id==rootID).and.(iter.eq.1)) print *, K-1, DABS(DYE)/INITYE, ct2a-ct1a

 10      end do                 ! End Jacobi diagonalization

         if (id==rootID) ct2s = MPI_WTIME()
         if ((id==rootID).and.(iter.eq.1)) print *
         if ((id==rootID).and.((iter.eq.1).or.(iter.eq.BLK))) print *, '   HOSB MPI Sum:',ITER, BLK, K-1, ct2s-ct1s
         !if (id==rootID) print *

         do I=1, MAK  
            RES(1:NDK,I,ITER) = HOSB(1:NDK,I,I,ITER)
         end do         

         do I = 1, MAK
            HOSB(1:NDK, i, i, iter) = 0.0D0
         end do
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         eff0 = abs(1.0D0 - ye/initye)
         if (iter==1) then
             eff(1) = eff0;     eff(2) = eff0
         else
             if (eff0 < eff(1) ) eff(1) = eff0
             if (eff0 > eff(2) ) eff(2) = eff0
         end if

         eff(3) = eff(3) + eff0
         eff(4) = eff(4) + sumDiag_MPI(COMM,NDK, MAK, HOSB(1,1,1,iter), IERR)
         eff(5) = eff(5) + k
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      end do MAIN_LOOP    

!      if(ID==rootID) then
!          print *
!          write(*, 20) eff(5)/BLK,  eff(3)/eff(5)
!          write(*, 30) eff(1), eff(2), eff(3)/BLK
!          print *,  '----------------------------------------------------------------------'
!      end if

20        format('  Aver. # of Sweeps :', F10.5, ',   Aver. OffDiag/Diag :',F10.6)
30        format('  Jacobi Effi: Min :', F10.5, ',   Max :', F10.5,',  Mean :',F10.5)
  end
!****************************************************************************

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                      SumOffDiag function:                                c
!c     Calculates the total error in the off-block-diagonal elements of H   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  double precision function sumOffDiag_MPI(COMM, NIN,NOUT,HOUT,IERR)
      implicit none
      include 'mpif.h'

      integer, intent(IN) :: COMM, NIN, NOUT             
      double precision, intent(IN) :: Hout(NIN,NOUT,NOUT)
      integer, intent(OUT) :: IERR

      integer :: I, J       
      double precision :: locSum

      locSum = 0.0D0
      do I = 1,NOUT  
         do J = 1,i-1
            locSum = locSum+dot_product(HOUT(1:NIN,J,I),HOUT(1:NIN,J,I)) 
         end do
      end do

      call MPI_ALLREDUCE(locSum, sumOffDiag_MPI,1,MPI_DOUBLE_PRECISION, &
             MPI_SUM, COMM, ierr)

      sumOffDiag_MPI = 2.0d0 * sumOffDiag_MPI

  end
!********************************************************************

  double precision function sumDiag_MPI(COMM, NIN,NOUT,HOUT,IERR)
      implicit none
      include 'mpif.h'

      integer, intent(IN) :: COMM, NIN, NOUT             
      double precision, intent(IN) :: Hout(NIN,NOUT,NOUT)
      integer, intent(OUT) :: IERR

      integer :: I
      double precision :: locSum

      locSum = 0.0D0
      do I = 1,NOUT  
         locSum = locSum+dot_product(HOUT(1:NIN,I,I),HOUT(1:NIN,I,I))
      end do

      call MPI_ALLREDUCE(locSum, sumDiag_MPI,1,MPI_DOUBLE_PRECISION, &
             MPI_SUM, COMM, ierr)    

  end
!********************************************************************
