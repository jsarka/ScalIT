!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                     c
!c     Subroutine to the abscissas and weights for     c
!c   Gauss-Laguerre function: See "Numerical Recipes"  c
!c                                                     c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Get the DVR points for Associated Laguerre Polynimial      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function AssLaguX0(alpha, N, x0)
   implicit none
   double precision, intent(IN)  :: alpha
   integer, intent(IN)  :: N
   double precision, intent(OUT) :: X0(N)

   double precision :: beta(N),vec, work
   integer :: i, ldz, info

   ldz = 1;  X0(1:N)=0.0D0
   DO I = 1, N
       x0(I) = 2.0D0*I -1.0D0 + alpha
       beta(I) = -DSQRT((alpha+i)*i)
   END DO

   CALL DSTEV('N', N, X0, beta, vec, ldz, work, info)

   AssLaguX0 = (info==0)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function AssLaguX0W0(alpha, N, x0, w0)
   implicit none
   double precision, intent(IN)  :: alpha
   integer, intent(IN)  :: N
   double precision, intent(OUT) :: X0(N),W0(N)

   double precision :: beta(N), vec(N,N), work(3*N)
   double precision :: gammln, nc0
   integer :: i, ldz, info

   ldz = N;  X0(1:N)=0.0D0
   DO I = 1, N
       x0(I) = 2.0D0*I -1.0D0 + alpha
       beta(I) = -DSQRT((alpha+i)*i)
   END DO

   CALL DSTEV('V', N, X0, beta, vec, ldz, work, info)

   AssLaguX0W0 = (info==0)

   if (AssLaguX0W0) then
       nc0=gammln(1.0+alpha)
       W0(1:N)=exp(nc0)*vec(1,1:N)*vec(1,1:N)
   end if
end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function AssLaguX0Vec(alpha, N, x0, Vmat)
   implicit none
   double precision, intent(IN)  :: alpha
   integer, intent(IN)  :: N
   double precision, intent(OUT) :: X0(N), Vmat(N,N)

   double precision :: beta(N), work(3*N)
   integer :: i, ldz, info

   ldz = N;  X0(1:N)=0.0D0
   DO I = 1, N
       x0(I) = 2.0D0*I -1.0D0 + alpha
       beta(I) = -DSQRT((alpha+i)*i)
   END DO

   CALL DSTEV('V', N, X0, beta, Vmat, ldz, work, info)

   AssLaguX0Vec = (info==0)

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine GauAssLaguX0W0(alpha, N,X0, W0)
   implicit none
   double precision, intent(IN)  :: alpha
   integer, intent(IN)  :: N
   double precision, intent(OUT) :: X0(N),W0(N)

   call gau_Laguerre(alpha,N,X0,W0)

end 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
