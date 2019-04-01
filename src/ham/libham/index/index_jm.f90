!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Index between k <-> (j,m)            c
!c  Here, j = 0, 1, ..., jmax                      c
!c        m = -j, -j+1, ..., -1, 0, 1, ..., j-1, j c 
!c                           for (j1m1)            c
!c        m = 0, 1, ..., j  for (j2m2)             c
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!c Total size of indices (k/(jm)) for j up to jmax   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function  getJm1Size(jmax)
    implicit none
    integer, intent(IN) :: jmax

    getJm1Size = (jmax+1)**2
end

integer function  getJm2Size(jmax)
    implicit none
    integer, intent(IN) :: jmax

    getJm2Size = (jmax+2)*(jmax+1)/2
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccc 
!c                (j, M) -> k                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function getJm1Pos(j, m)
    implicit none
    integer, intent(IN)  ::  j, m
    
    getJm1Pos = j**2 + j + m + 1
end

integer function getJm2Pos(j, m)
    implicit none
    integer, intent(IN)  ::  j, m
    
    getJm2Pos = (j+1)*j/2 + m + 1
end

!cccccccccccccccccccccccccccccccccccccccccccccccccc
!c                 k -> (j, M)                    c
!cccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getJm1Index(k, j, m)
   implicit none
   integer, intent(IN)  :: k
   integer, intent(OUT) :: j, m

   integer :: j0, j1, num0, num1, ind

   j0 = 0; j1=j0+1

   do        
       num0 = j0**2;    num1 = j1**2
       if ((k>num0) .AND. (k<=num1)) then
           j = j0; exit
       else
           j0=j1; j1=j1+1
       end if
    end do
 
    m = k - j**2 - j -1
end 

subroutine getJm2Index(k, j, m)
   implicit none
   integer, intent(IN)  :: k
   integer, intent(OUT) :: j, m

   integer :: j0, j1, num0, num1, ind

   j0 = 0; j1=j0+1

   do
       num0 = (j0+1)*j0/2;    num1 = (j1+1)*j1/2
       if ((k>num0) .AND. (k<=num1)) then
           j = j0; exit
       else
           j0=j1; j1=j1+1
       end if
    end do
 
    m = k - (j+1)*j/2 - 1
end 


!cccccccccccccccccccccccccccccccccccccccccccccccccc
!c           k -> (N(k), L(k), M(k))              c
!cccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getJm1Indices(kmax, jm)
   implicit none
   integer, intent(IN)  :: kmax
   integer, intent(OUT) :: jm(2, kmax)

   integer :: jmax, j, M, k

   jmax = sqrt(dble(kmax)) + 3

   k = 0
   do j = 0, jmax
      do M = -j, j
         k = k + 1
         if (k > kmax) then
            return
         else
            jm(1,k)=j; jm(2,k)=m
         end if 
      end do
   end do

end

subroutine getJm2Indices(kmax, jm)
   implicit none
   integer, intent(IN)  :: kmax
   integer, intent(OUT) :: jm(2, kmax)

   integer :: jmax, j, M, k

   jmax = sqrt(dble(kmax*2)) + 3

   k = 0
   do j = 0, jmax
      do M = 0, j
         k = k + 1
         if (k > kmax) then
            return
         else
            jm(1,k)=j; jm(2,k)=m
         end if    
      end do
   end do

end

!*********************************************************
