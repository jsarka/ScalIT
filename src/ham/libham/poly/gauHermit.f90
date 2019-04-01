!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Subroutine to the abscissas and weights for       c
!c   Gauss-Hermite function: See "Numerical Recipes" 	 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine HermiteNode(N, Xout)
   implicit none
   integer, intent(IN)  :: N          
   double precision, intent(OUT) :: Xout(N)   

   double precision, dimension(N)  :: beta
   double precision :: vec, work
   integer :: i, ldz, info

   ldz = 1

   DO I = 1, N
       beta(I) =  DSQRT(0.5D0*I)
   END DO

   CALL DSTEV('N', N, XOUT, beta, vec, ldz, work, INFO)

end

subroutine hermite_alpha(alpha, N, Xout)
    !*******************************************************
   implicit none
   double precision, intent(IN) :: alpha
   integer, intent(IN)          :: N          
   double precision, intent(OUT), dimension(N) :: Xout   

   double precision, dimension(N)  :: beta
   double precision :: vec, work
   integer :: i, ldz, info

   ldz = 1

   DO I = 1, N
       beta(I) =  DSQRT(0.5D0*I)
   END DO

   CALL DSTEV('N', N, XOUT, beta, vec, ldz, work, INFO)

   Xout(1:N) = XOUT(1:)/alpha    ! alpha=sqrt(m*freq/h)
end


!********************************************************************
SUBROUTINE GAU_HERMITE(n, x, w)
   implicit none
   double precision, parameter  :: EPSI = 1.0D-14
   ! double precision, parameter  :: PIE  = 3.14159265358979324D0
   double precision, parameter  :: PIM4 = .7511255444649425D0   ! 1/pi^(1/4)
   integer, parameter :: MAXIT = 100 

   integer, intent(in)          :: n       !number of abscissas and weights
   double precision, intent(out), dimension(n) :: x  !real abscissas 
   double precision, intent(out), dimension(n) :: w  !real weights
   
   integer :: i,its,j,m
   double precision :: p1,p2,p3,pp,z,z1

    m=(n+1)/2
    do i=1,m
        if(i == 1)then
          z = dsqrt(1.0D0*(2*n+1))-1.85575D0/(2*n+1)**(0.16667D0)
        else if(i == 2)then
          z=z-1.14*n**.426/z
        else if (i == 3)then
          z=1.86*z - .86*x(1)
        else if (i == 4)then
          z=1.91*z-.91*x(2)
        else
          z=2.*z-x(i-2)
        endif

        do its=1,MAXIT
          p1=PIM4
          p2=0.d0

          do j=1,n
            p3=p2
            p2=p1
            p1=z*dsqrt(2.d0/j)*p2-dsqrt(dble(j-1)/dble(j))*p3
          end do

          pp = dsqrt(2.d0*n)*p2
          z1=z
          z=z1-p1/pp

          if(abs(z-z1) <= EPSI) exit

        end do

        x(i)=z
        x(n+1-i)=-z
        w(i)=2.d0/(pp*pp)
        w(n+1-i)=w(i)

   end do
      
END SUBROUTINE 
