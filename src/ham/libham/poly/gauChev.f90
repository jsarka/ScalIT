!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Subroutine to the abscissas and weights for       c
!c   Gauss-Chebychev function: See "Numerical Recipes" 	 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine chevNode(N, Xout)  
   implicit none
   integer, intent(IN)  :: N          
   double precision, intent(OUT) :: Xout(N)   

   double precision, dimension(N)  :: beta
   double precision :: vec, work
   integer :: i, ldz, info

   ldz = 1

   do I = 1, N
       beta(I) =  0.5D0
   end do

   CALL DSTEV('N', N, XOUT, beta, vec, ldz, work, INFO)

end


!********************************************************************
SUBROUTINE gau_chev(n, x, w)
   implicit none
   integer, intent(in)          :: n       !number of abscissas and weights
   double precision, intent(out) :: x(N), w(N)

   double precision, parameter  :: PIE  = 3.14159265358979324D0   

   integer :: i

   DO I = 1, N
      x(I) = DCOS(PIE*(I-0.5D0)/N)
   END DO
   
   W(1:N) = PIE/N
      
END 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
