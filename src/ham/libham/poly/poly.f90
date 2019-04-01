!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Polynomials. Calculate the polynomial:      c
!c    Y=C(1)*X^N+C(2)*X^(N-1)+...+C(N)*X+C(N+1)       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine PolyVal(N, coeff, M, X, Y)
   implicit none
   integer, intent(IN) :: N, M
   double precision, intent(IN) :: coeff(N+1), X(M)
   double precision, intent(OUT):: Y(M)

   integer :: i, j

   Y(1:M) = coeff(1)*X(1:M)
   do i = 2, N
      Y(1:M) = (Y(1:M)+coeff(i))*X(1:M)
   end do
   Y(1:M) = Y(1:M) + coeff(N+1)
end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Calculate the step tanh(x) function              c
!c  Y = coeff(1)* tanh(coeff(2)*(x-coeff(3))) + coeff(4)  c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine myStep(coeff, M, X, Y)
   implicit none
   integer, intent(IN) :: M
   double precision, intent(IN) :: coeff(4), X(M)
   double precision, intent(OUT):: Y(M)

   Y(1:M) = coeff(1)*tanh(coeff(2)*(X(1:M)-coeff(3)))+coeff(4);

end 

!cccccccccccccccccccccccccccccccccccccccccccccc
!c     Calculate the mid-filter function      c
!cccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine myMidFilter(coeff, M, X, Y)
   implicit none
   integer, intent(IN) :: M
   double precision, intent(IN) :: coeff(6), X(M)
   double precision, intent(OUT):: Y(M)

   double precision :: X1(M), X2(M)

   X1(1:M) = coeff(2)*(X(1:M)-coeff(3))
   X2(1:M) = coeff(4)*(X(1:M)-coeff(5))
   Y(1:M) = coeff(1)*tanh(X1(1:M))*tanh(X2(1:M))+coeff(6);

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Calculate the Gaussian-Type oscillator function     c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine myGOSC(coeff, M, X, Y)
   implicit none
   integer, intent(IN) :: M
   double precision, intent(IN) :: coeff(5), X(M)
   double precision, intent(OUT):: Y(M)

   Y(1:M) =  coeff(1) * cos(coeff(2)*(X(1:M)-coeff(3)))     &
           * exp(-coeff(4)*(X(1:M)-coeff(3))**2) + coeff(5)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

