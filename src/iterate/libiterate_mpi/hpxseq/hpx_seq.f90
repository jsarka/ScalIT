!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c          H1_MULT_X subroutine, Real version                c
!c     Y(nin,nout)=H(nout,nout)@I(nin,nin)*X(nin,nout)        c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine H1X_XYZ_Seq(NIN, NOUT, BLK, H1, X, Y)
      implicit none
      integer,intent(IN) :: NIN, NOUT, BLK
      double precision,intent(IN)  :: H1(NOUT,NOUT)  
      double precision,intent(IN)  :: X(NIN,NOUT,BLK)  
      double precision,intent(OUT) :: Y(NIN,NOUT,BLK)

      integer :: i, j, k

      do k = 1, BLK

         do i = 1, NOUT

            Y(1:NIN, i, k) = H1(i, 1) * X(1:NIN, 1, k)

            do j=2,NOUT

               Y(1:NIN,i,k) = Y(1:NIN,i,k)+H1(i,j)*X(1:NIN,j,k)  

            end do

         end do

      end do

 end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine H1TX_XYZ_Seq(NIN, NOUT, BLK, H1, X, Y)
      implicit none
      integer,intent(IN) :: NIN, NOUT, BLK
      double precision,intent(IN)  :: H1(NOUT,NOUT)  
      double precision,intent(IN)  :: X(NIN,NOUT,BLK)  
      double precision,intent(OUT) :: Y(NIN,NOUT,BLK)

      integer :: i, j, k

      do k = 1, BLK

         do i = 1, NOUT

            Y(1:NIN, i, k) = H1(1,i) * X(1:NIN, 1, k)

            do j=2, NOUT

               Y(1:NIN,i,k) = Y(1:NIN,i,k)+H1(j,i)*X(1:NIN,j,k)  

            end do

         end do

      end do
 end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!**************************************************************

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine H2X_DEP_Seq(NIN, NOUT, BLK, H1, Depx, X, Y)
      implicit none
      integer,intent(IN) :: NIN, NOUT, BLK
      double precision,intent(IN)  :: H1(NOUT,NOUT)       
      double precision,intent(IN)  :: DepX(NIN)   
      double precision,intent(IN)  :: X(NIN,NOUT,BLK)  
      double precision,intent(OUT) :: Y(NIN,NOUT,BLK)

      integer :: i, j, k

      do k = 1, BLK

         do i = 1, NOUT

            Y(1:NIN, i, k) = H1(i,1)*X(1:NIN,1,k)*DepX(1:NIN)

            do j=2,NOUT

               Y(1:NIN,i,k) = Y(1:NIN,i,k)+H1(i,j)*X(1:NIN,j,k)*DepX(1:NIN)

            end do

         end do

      end do
 end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine H2TX_DEP_Seq(NIN, NOUT, BLK, H1, Depx, X, Y)
      implicit none
      integer,intent(IN) :: NIN, NOUT, BLK
      double precision,intent(IN)  :: H1(NOUT,NOUT)       
      double precision,intent(IN)  :: DepX(NIN)   
      double precision,intent(IN)  :: X(NIN,NOUT,BLK)  
      double precision,intent(OUT) :: Y(NIN,NOUT,BLK)

      integer :: i, j, k

      do k = 1, BLK

         do i = 1, NOUT

            Y(1:NIN, i, k) = H1(1,i)*X(1:NIN,1,k)*DepX(1:NIN)

            do j=2,NOUT

               Y(1:NIN,i,k) = Y(1:NIN,i,k)+H1(j,i)*X(1:NIN,j,k)*DepX(1:NIN)

            end do

         end do

      end do

 end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!**************************************************************

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine H3X_Out_Seq(NIN, NOUT, BLK, outH, X, Y)
      implicit none
      integer,intent(IN) :: NIN, NOUT, BLK
      double precision,intent(IN)  :: outH(NIN,NOUT,NOUT)  
      double precision,intent(IN)  :: X(NIN,NOUT,BLK)  
      double precision,intent(OUT) :: Y(NIN,NOUT,BLK)

      integer :: i, j, k

      do k = 1, BLK

         do i = 1, NOUT

            Y(1:NIN, i, k) = outH(1:NIN, i, 1) * X(1:NIN, 1, k)

            do j=2,NOUT

               Y(1:NIN,i,k) = Y(1:NIN,i,k)+outH(1:NIN,i,j)*X(1:NIN,j,k)  

            end do

         end do

      end do

 end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine H3TX_Out_Seq(NIN, NOUT, BLK, outH, X, Y)
      implicit none
      integer,intent(IN) :: NIN, NOUT, BLK
      double precision,intent(IN)  :: outH(NIN,NOUT,NOUT)  
      double precision,intent(IN)  :: X(NIN,NOUT,BLK)  
      double precision,intent(OUT) :: Y(NIN,NOUT,BLK)

      integer :: i, j, k

      do k = 1, BLK

         do i = 1, NOUT

            Y(1:NIN, i, k) = outH(1:NIN,1,i) * X(1:NIN, 1, k)

            do j=2, NOUT

               Y(1:NIN,i,k) = Y(1:NIN,i,k)+outH(1:NIN,j,i)*X(1:NIN,j,k)  

            end do

         end do

      end do

 end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      V_MULT_X: Do the real V*X job at specific level       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine VX_Seq(NIN, NOUT, BLK, VOUT, X)
      implicit none      
      integer, intent(IN) :: NIN, NOUT, BLK
      double precision,intent(IN)    :: VOUT(NOUT,NOUT,BLK) 
      double precision,intent(INOUT) :: X(NIN,NOUT, BLK)

      integer  :: i, j, k
      double precision :: TMP(NIN,NOUT) 

      do k = 1, BLK

         do i = 1, NOUT

            tmp(1:NIN, i) = Vout(i, 1, k) * X(1:NIN, 1, k)

            do j=2,NOUT

               tmp(1:NIN,i) = tmp(1:NIN,i)+Vout(i,j,k)*X(1:NIN,j,k)  

            end do

         end do

         X(1:NIN, 1:NOUT, k) = TMP(1:NIN, 1:NOUT)

      end do

 end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     VT_MULT_X: Do the real V^T*X job at specific level       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine VTX_Seq(NIN, NOUT, BLK, VOUT, X)
      implicit none      
      integer, intent(IN) :: NIN, NOUT, BLK
      double precision, intent(IN) :: Vout(NOUT,NOUT,BLK)
      double precision, intent(INOUT) :: X(NIN,NOUT, BLK)

      integer ::  i, j, k
      double precision :: TMP(NIN,NOUT)     

      do k = 1, BLK

         do i = 1, NOUT

            tmp(1:NIN, i) = Vout(1, i, k) * X(1:NIN, 1, k)

            do j=2,NOUT

               tmp(1:NIN,i) = tmp(1:NIN,i)+Vout(j,i,k)*X(1:NIN,j,k)  

            end do

         end do

         X(1:NIN, 1:NOUT, k) = TMP(1:NIN, 1:NOUT)

      end do

 end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
