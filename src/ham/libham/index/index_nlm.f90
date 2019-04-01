!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Index between k <-> (n,l,m)          c
!c  Here, N = 1, 2, ..., Nmax; L = 0, 1, ..., N-1  c
!c        M = -L, -L+1, ..., -1, 0, 1, ..., L-1, L c 
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!c Total size of indices (k/(nlm)) for N up to Nmax  c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function  getNLMSize(nmax)
    implicit none
    integer, intent(IN) :: nmax

    getNLMSize = (nmax**3 + (3*nmax**2 + nmax)/2) /3
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccc 
!c              (N, L, M) -> k                     c
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function getNLMPos(N,L,M)
    implicit none
    integer, intent(IN)  ::  N, L, M
    
    getNLMPos=((N-1)**3+(3*(N-1)**2+(N-1))/2)/3 + (L+1)**2 - L + M;
end

!cccccccccccccccccccccccccccccccccccccccccccccccccc
!c                k -> (N, L, M)                  c
!cccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getNLMIndex(k, N, L, M)
   implicit none
   integer, intent(IN)  :: k
   integer, intent(OUT) :: N, L, M

   integer :: n0, n1, num0, num1, ind

   n0 = 0; n1=n0+1

   do        
       num0 = (n0**3+(3*n0**2+n0)/2)/3
       num1 = (n1**3+(3*n1**2+n1)/2)/3
       if ((k>num0) .AND. (k<=num1)) then
           N = n1; exit
       else
           n0=n1; n1=n1+1
       end if
    end do
 
    n1 = N-1
    ind = k-(n1**3 + (3*n1**2+n1)/2)/3;
    L = floor(sqrt(ind-1.0D0))
    M = ind - L**2 - l - 1
end 

!cccccccccccccccccccccccccccccccccccccccccccccccccc
!c           k -> (N(k), L(k), M(k))              c
!cccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getNLMIndices(kmax, NLM)
   implicit none
   integer, intent(IN)  :: kmax
   integer, intent(OUT) :: NLM(3, kmax)

   integer :: nmax, N, L, M, k

   nmax = (3*kmax)**(1.0D0/3.0D0) + 2  ! 1 

   k = 0
   do N = 1, nmax
      do L = 0, N-1
         do M = -L, L
            k = k + 1
            if (k > kmax) then
                return
            else
                NLM(1,k)=N; NLM(2,k)=L;NLM(3,k)=M
            end if 
         end do
      end do
   end do

end

!*********************************************************
