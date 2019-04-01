!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c          Index for Tetra-atomic molecules: 1D<->(j1j2jK)                 c
!c          !! MODIFIED BY BILL POIRIER ON 4/28/14 !!                       c
!c                                                                          c
!c This is for A2B2 molecules with G4 Aug permutation parity for (j1j2):    c
!c Aug: j1=odd, j2=even                                                     c
!c                                                                          c
!c parity is TRUE parity, including (-1)^J factor, what Zhang calls "p".    c
!c          in tetra-atomic case, parity adaptation is thus:                c
!c          |K> + (-1)^(p+j1+j2+j+J)|-K>                                    c
!c (j,k) are combined angular momentum functions from (j1,m1) and (j2,m2)   c
!c JTol is big J (total angular momentum, J)                                c  
!c jmax does additional basis truncation (if jmax < j1max+j2max)            c
!c jksize is combined angular basis size (from 4 indices)                   c
!c mmsize is total number of Clebsch-Gordon coefficients for all (j,k)'s    c 
!c i0 is the composite index for the angular functions, (j1,j2,j,k).        c
!c jkind stores the (j1,j2,j,k) values associated with each i0.             c
!c msize is the number of CG coeffs for each i0 angular function.           c
!c mindex stores all of the individual m values for all CG coeffs.          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine get4Size(parity, JTol, j1max,j2max,jmax, jkSize, mmSize)
   implicit none
   logical, intent(IN)  :: parity
   integer, intent(IN)  :: JTol, j1max, j2max, jmax
   integer, intent(OUT) :: jkSize, mmSize

   integer :: j1, j2, j, k, m, km, kmin, jsum
   integer :: j0min, j0max, k0max

   jkSize = 0; mmSize = 0
   do j1 = 1, j1max, 2
      do j2 = 0, j2max, 2
         j0min = ABS(j1-j2)
         j0max = min(j1+j2, jmax)
         do j = j0min, j0max
              k0max = min(j,Jtol)
              jsum = j1 + j2 + j + JTol
              if (.NOT. parity) jsum = jsum + 1
              if (jsum/2*2 == jsum) then  ! include k=0
                 kmin = 0
              else
                 kmin = 1
              end if
              jkSize = jkSize + k0max + 1 - kmin
              do k=kmin, k0max
                  do m=-j1, j1
                     km = ABS(k-m)
                     if (.NOT.(km>j2)) mmSize = mmSize + 1 
                  end do
              end do
         end do
      end do
   end do
end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine get4Index(parity, JTol, j1max,j2max,jmax, N, jkind, mSize, N0, mIndex)
   implicit none
   logical, intent(IN)  :: parity
   integer, intent(IN)  :: JTol, j1max, j2max, jmax,  N, N0
   integer, intent(OUT) :: jkInd(4,N), mSize(N), mIndex(N0)

   integer :: j1, j2, j, k, m, km, i0, m0, kmin, jsum
   integer :: j0min, j0max, k0max

   i0 = 0; m0=0
   do j1 = 1, j1max, 2
      do j2 = 0, j2max, 2
         j0min = ABS(j1-j2)
         j0max = min(j1+j2, jmax)
         DO j = j0min, j0max
              k0max = min(j, JTol)
              jsum = j1 + j2 + j + JTol
              if (.NOT. parity) jsum = jsum + 1
              if (jsum/2*2 == jsum) then  ! include k=0
                 kmin = 0
              else
                 kmin = 1
              end if
              DO k = kmin, k0max 
                 i0 = i0 + 1
                 if (i0>N) return  ! angular basis larger than some allowed value N 
                 jkInd(1,i0) = j1; jkind(2,i0)=j2
                 jkind(3,i0) = j;  jkind(4,i0)=k
                 mSize(i0) = 0;
                 DO m = -j1, j1
                     km = ABS(K-m)
                     if (.NOT. (km>j2)) then
                        m0 = m0 + 1
                        if (m0>N0)  return  ! number of CG coeffs larger than some allowed N0
                        mSize(i0)= mSize(i0)+1 
                        mIndex(m0) = m    
                     end if
                 END DO  ! m
              END DO     ! k
         END DO          ! j
      end do  ! j2
   end do     ! j1
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

