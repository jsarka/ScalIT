!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Subroutines to initialize HOSB and diagonize HOSB                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C     HOSBDIAG ROUTINES                                                       C
!C     The contract version to Perform Block-Jacobi Diagonization              c 
!C     After Block-Jacobi transformation,the transformation matrices           C
!C     are stored in VOSB,and HOSB reveals the results after the Block-Jacobi  C
!C     transformation,diagonal submatrices of HOSB are stored back to EIGVAL   C
!C          HOSB:  The initial H matrix at LEV=k,                              C
!C          VOSB:  The address to store transformation matrix.                 C
!C          EIGVALL:   The diagonal submatrices.                               C
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function BJDIAG( BJMAX, BJTOL, NIN, NOUT, HOSB, VOSB, EIGVAL)
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: BJMAX, NIN, NOUT
      DOUBLE PRECISION, INTENT(IN) :: BJTOL 

      DOUBLE PRECISION,INTENT(INOUT):: Hosb(NIN, NOUT, NOUT)
      DOUBLE PRECISION,INTENT(OUT)  :: VOSB(NOUT, NOUT)
      DOUBLE PRECISION,INTENT(OUT):: EIGVAL(NIN, NOUT)    

!ccccccccccccccccccccccccccccccc Local variables
      DOUBLE PRECISION :: ye, initye,dye    
      DOUBLE PRECISION :: sumOffDiag, PHI, MJACOBI   
      INTEGER :: i, j, k       

      initye = SumOffDiag(NIN, NOUT,HOSB(1,1,1))         
         
      ye = initye;   dye = initye ;       k = 1

      ! initial block NOUT x NOUT of V to I
      CALL VINIT(NOUT, VOSB)

      DO 10 WHILE((DABS(dye) > initye * BJTOL).AND.(k <= BJMAX))
          DO i = 1,(NOUT-1)   
               DO j = 1,(NOUT-I) 
                  PHI = MJACOBI(j, j+i,NOUT,NIN, HOSB(1,1,1))
                  CALL VUPDATE(j, j+i, PHI, NOUT, VOSB)
               END DO
          END DO              ! One Jacobi cycle                      

          dye = SumOffDiag(NIN,NOUT,HOSB(1,1,1)) - ye
          ye = ye + dye
            
          k = k + 1            

 10   END DO                 ! End Jacobi diagonalization

!********************************************************
      BJDIAG = ABS(1.0D0 - ye/initye)

      DO i=1, NOUT            
         EIGVAL(1:NIN,i) = HOSB(1:NIN, i, i)
      END DO                    

END 
!******************************************************************************

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                      YERROR function:                                    c
!c     Calculates the total error in the off-block-diagonal elements of H   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 DOUBLE PRECISION FUNCTION SumOffDiag(NIN,NOUT,HOUT)
      IMPLICIT NONE
      INTEGER,INTENT(IN)   ::  NIN, NOUT    
      DOUBLE PRECISION,INTENT(IN) :: HOUT(NIN,NOUT,NOUT)    
      
      INTEGER  ::  i, j            
   
      SumOffDiag = 0.0
      DO i = 1,(NOUT-1)        
         DO j = 1,(NOUT-I)
            SumOffDiag = SumOffDiag + DOT_PRODUCT(HOUT(1:NIN,j,j+i), HOUT(1:NIN,j,j+i))
         END DO
      END DO

      SumOffDiag = 2.0d0 * SumOffDiag

END
!***************************************************************************

