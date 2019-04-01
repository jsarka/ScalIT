

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Summation of the off-diagonal blocks                   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 double precision function sumOffDiag(NIN,NOUT,HOUT)
      implicit none
      integer,intent(IN)           :: NIN, NOUT   
      double precision, intent(IN) :: Hout(NIN,NOUT,NOUT)    
      
      integer  ::  i, j            
   
      sumOffDiag = 0.0
      do i = 1,NOUT     
         do j = 1, i-1
           sumOffDiag = sumOffDiag + dot_product(HOUT(1:NIN,j,i),HOUT(1:NIN,j,i))
         end do
      end do
  
      sumOffDiag = 2.0D0*sumOffDiag

 end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Summation of the diagonal blocks                       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 double precision function sumDiag(NIN,NOUT,HOUT)
      implicit none
      integer,intent(IN)           :: NIN, NOUT   
      double precision, intent(IN) :: Hout(NIN,NOUT,NOUT)    
      
      integer  ::  i
   
      sumDiag = 0.0
      do i = 1,NOUT     
          sumDiag = sumDiag + dot_product(HOUT(1:NIN,i,i), HOUT(1:NIN,i,i))
      end do

 end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Summation of the off-diagonal blocks                   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 double precision function sumOffDiagCX(NIN,NOUT,HOUT)
      implicit none
      integer,intent(IN)           :: NIN, NOUT   
      double precision, intent(IN) :: Hout(NIN,NOUT,NOUT)    
      
      integer  ::  i, j            
   
      sumOffDiagCX = 0.0
      do i = 1,NOUT     
         do j = 1, i-1
           sumOffDiagCX = sumOffDiagCX + ABS(dot_product(HOUT(1:NIN,j,i),   &
                                         HOUT(1:NIN,j,i)))
         end do
      end do
  
      sumOffDiagCX = 2.0D0*sumOffDiagCX
 end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Summation of the diagonal blocks                       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 double precision function sumDiagCX(NIN,NOUT,HOUT)
      implicit none
      integer,intent(IN)           :: NIN, NOUT   
      double precision, intent(IN) :: Hout(NIN,NOUT,NOUT)    
      
      integer  ::  i
   
      sumDiagCX = 0.0
      do i = 1,NOUT     
          sumDiagCX = sumDiagCX + ABS(dot_product(HOUT(1:NIN,i,i),       &
                                    HOUT(1:NIN,i,i)))
      end do

 end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

