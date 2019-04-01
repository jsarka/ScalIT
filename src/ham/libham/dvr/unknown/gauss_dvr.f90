!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   DVR functions obtained from Gaussian basis set  c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine Gauss_XS(N, a, b, alpha, S1, X1)
   implicit none
   integer, intent(IN) :: N
   double precision, intent(IN) :: a, b, alpha
   double precision, intent(OUT):: S1(N,N), X1(N,N)

   integer :: i, j, ij
   double precision :: dx, c0, x0(N),ex0(N)

   dx=(b-a)/(N-1); c0=0.5D0*alpha*dx*dx;
   do i = 1, N
      x0(i)=a+(i-1)*dx
      ex0(i)=exp(-c0*i*i);
   end do

   do i = 1, N
      do j = 1,i-1 
         ij=i-j;    S1(i,j)=ex0(ij)
         X1(i,j)=0.5D0*(x0(i)+x0(j))*S1(i,j)
         S1(j,i)=S1(i,j)
         X1(j,i)=X1(i,j)
      end do
      S1(i,i)=1.0;
      X1(i,i)=x0(i)
   end do

end 

subroutine GaussCx_XS(N, a, b, alpha, S1, X1)
   implicit none
   integer, intent(IN) :: N
   double precision, intent(IN) :: a, b, alpha
   double complex, intent(OUT):: S1(N,N), X1(N,N)

   integer :: i, j, ij
   double precision :: dx, c0, x0(N),ex0(N)

   dx=(b-a)/(N-1); c0=0.5D0*alpha*dx*dx;
   do i = 1, N
      x0(i)=a+(i-1)*dx
      ex0(i)=exp(-c0*i*i);
   end do

   do i = 1, N
      do j = 1,i-1 
         ij=i-j;    S1(i,j)=ex0(ij)
         X1(i,j)=0.5D0*(x0(i)+x0(j))*S1(i,j)
         S1(j,i)=Conjg(S1(i,j))
         X1(j,i)=Conjg(X1(i,j))
      end do
      S1(i,i)=1.0;
      X1(i,i)=x0(i)
   end do

end 
   

