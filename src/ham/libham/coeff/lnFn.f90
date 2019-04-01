!cccccccccccccccccccccccccccccccccccccccccccccc
!c         store ln(n!) in an array           c
!cccccccccccccccccccccccccccccccccccccccccccccc
subroutine lnFn(N, fn)          ! store 0!, 1!, ..., (N-1)!, fn(i) = ln(i-1)!
   implicit none
   integer, intent(IN) :: N
   double precision, intent(OUT) :: fn(N)

   integer :: i  
 
   fn(1:2) = 0.0D0

   do i=3, N
      fn(i) = fn(i-1) + LOG(DBLE(i-1))
   end do

end

subroutine lnFn1(N, fn)          ! store 1!, ..., (N)!,fn(i)=ln(i!)
   implicit none
   integer, intent(IN) :: N
   double precision, intent(OUT) :: fn(N)

   integer :: i  
 
   fn(1:2) = 0.0D0

   do i=2, N
      fn(i) = fn(i-1) + LOG(DBLE(i))
   end do

end

  

