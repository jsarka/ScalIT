!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   The transform between 1D and multi-dimension indices    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                           c
!c  (ind1(N))->get1DIndex, The dimensionality is ind0(1:N)   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function get1DLength(N, ind0)
   implicit none
   integer, intent(IN) :: N, ind0(N)

   integer :: i

   get1DLength = 1
   do i = 1, N
      get1DLength = get1DLength*ind0(i)
   end do

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function get1DIndex(N, ind0, ind1)
   implicit none
   integer, intent(IN) :: N, ind0(N), ind1(N)

   integer :: nDim(N), i

   nDim(1)=1; get1DIndex=ind1(1)
   do i = 2, N
      nDim(i) = nDim(i-1)*ind0(i-1)
      get1DIndex = get1DIndex+(ind1(i)-1)*nDim(i)
   end do

   print *, 'leave get1Dindex:', get1DIndex

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    (ind->ind1(N))    The dimensionality is ind0(1:N)      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getMDIndex(N, ind0, ind, ind1)
   implicit none
   integer, intent(IN)  :: N, ind0(N), ind
   integer, intent(OUT) :: ind1(N)

   integer :: nDim(N), i, dex

   nDim(1)=1
   do i=1, N-1
      nDim(i+1)=nDIm(i)*ind0(i)
   end do

!  print *, 'nDIm:', nDim(1:N)

   dex = ind-1
   do i=N,2, -1
      ind1(i)= dex/nDim(i)
      dex = dex - ind1(i)*nDim(i)
   end do

   ind1(1)=dex
   ind1(1:N) = ind1(1:N)+1
 !  ind1(1) = dex

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
