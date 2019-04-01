!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Subroutine to calculate sin(m*a) and cos(m*a)    c
!c               for the given angle a.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      cos(n*a), sin(n*a): n=1, 2, ..., Nmax           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine Trigo1(Nmax, N, theta, cnth, snth)
   implicit none
   integer, intent(IN) :: Nmax, N
   double precision, intent(IN) :: theta(N)
   double precision, intent(OUT):: cnth(N, Nmax), snth(N, Nmax)

   double precision :: cth(N), sth(N)

   cth(1:N) = cos(theta(1:N)); sth = sin(theta(1:N))

   call Trigo3(Nmax, N, cth, sth, cnTh, snTh)

end

subroutine Trigo2(Nmax, N, cth, cnth, snth)
   implicit none
   integer, intent(IN) :: NMax, N
   double precision, intent(IN) :: cth(N)
   double precision, intent(OUT):: cnth(N, NMax), snth(N, Nmax)

   double precision :: sth(N)

   sth(1:N) = sqrt(1.0D0 - cth(1:N)**2)

   call Trigo3(Nmax, N, cth, sth, cnTh, snTh)

end

subroutine Trigo3(Nmax, N, cth, sth, cnth, snth)
   implicit none
   integer, intent(IN) :: NMax, N
   double precision, intent(IN) :: cth(N), sth(N)
   double precision, intent(OUT):: cnth(N, NMax), snth(N, Nmax)

   integer :: i

   cnth(1:N, 1) = cth(1:N); snth(1:N, 1) = sth(1:N)
   cnth(1:N, 2) = 1.0D0 - 2.0D0*sth(1:N)*sth(1:N)
   snth(1:N, 2) = 2.0D0 * sth(1:N) * cth(1:N) 

   do i = 3, Nmax
       cnth(1:N, i) = 2.0D0*cnth(1:N, (i-1))*cth(1:N) - cnth(1:N, (i-2))
       snth(1:N, i) = 2.0D0*snth(1:N, (i-1))*cth(1:N) - snth(1:N, (i-2))
   end do

end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
