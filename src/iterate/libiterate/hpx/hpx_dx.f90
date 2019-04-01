!cccccccccccccccccccccccccccccccccccccccccccccccccc
!c       H1_MULT_X subroutine, Complex version    c
!c          Calculate Y = H(level) * X            c
!cccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine H1X_XYZ_DX(NIN, NOUT, BLK, H1, X, Y)
      implicit none
      integer,intent(IN) :: NIN, NOUT, BLK
      double precision,intent(IN) :: H1(NOUT,NOUT)  
      double complex,intent(IN)   :: X(NIN,NOUT,BLK)  
      double complex,intent(OUT)  :: Y(NIN,NOUT,BLK)  

      integer :: i

      do i = 1, BLK
         call H1X_DX(NIN, NOUT, H1, X(1,1,i), Y(1,1,i))
      end do

  end
!******************************************************************
  subroutine H1X_DEP_DX(NIN, NOUT, BLK, H1, CORX, X, Y)
      implicit none
      integer,intent(IN) :: NIN, NOUT, BLK
      double precision,intent(IN) :: H1(NOUT,NOUT)    
      double precision,intent(IN) :: CORX(NIN)     
      double complex,intent(IN)   :: X(NIN,NOUT,BLK)   
      double complex,intent(OUT)  :: Y(NIN,NOUT,BLK)

      integer :: i

      do i = 1, BLK
         call H2X_DX(NIN, NOUT, H1, CORX, X(1,1,i), Y(1,1,i))
      end do

  end
 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      V_MULT_X: Do the real V*X job at specific level       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine VX_DX(NIN, NOUT, BLK, VOUT, X)
      implicit none      
      integer, intent(IN) :: NIN, NOUT, BLK
      double precision,intent(IN)  :: Vout(NOUT,NOUT,BLK)
      double complex,intent(INOUT) :: X(NIN,NOUT, BLK) 

      integer  :: i 
      double complex :: TMP(NIN,NOUT) 

      do i = 1, BLK
         call H1X_DX(NIN, NOUT, VOUT(1,1,i), X(1,1,i), TMP)
         X(1:NIN, 1:NOUT, I) = TMP(1:NIN, 1:NOUT)
      end do

  end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     VT_MULT_X: Do the real V^T*X job at specific level       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine VTX_DX(NIN, NOUT, BLK, VOUT, X)
      implicit none      
      integer, intent(IN) :: NIN, NOUT, BLK
      double precision,intent(IN)  :: Vout(NOUT,NOUT,BLK)
      double complex,intent(INOUT) :: X(NIN,NOUT, BLK)

      integer ::  i
      double complex :: tmp(NIN,NOUT) 

      do i = 1, BLK
         call H1TX_DX(NIN, NOUT, VOUT(1,1,i), X(1,1,i), TMP)
         X(1:NIN, 1:NOUT, I) = TMP(1:NIN, 1:NOUT)
      end do

  end
!*****************************************************************

