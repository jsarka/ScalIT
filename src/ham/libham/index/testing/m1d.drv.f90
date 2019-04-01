!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   The grouping of nodes for each block           c
!c     This only works when node > blk              c
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program test_NodeLen
  integer, parameter :: NMAX=20
  integer :: ndim(NMAX), N, nind(NMAX), len, ind 
  integer :: get1DLength, get1DIndex

  print *, 'Input Number of dimension: N'
  read (*,*) N
  
  print *, ' Input sizes for each dimension:'
  read (*,*) ndim(1:N)


  print *, ' Dimension:', nDim(1:N)

  len = get1DLength(N, ndim)
  print *, ' Total Size:', len

  do i = 1, len
     call getMDIndex(N, nDim, i, nInd)
     print *, ' 1D Index:', i, ' -> MD Indices:', nInd(1:N)
    
     ind = get1DIndex(N, nDim, nind)   
     print *, ' MD Indices: ', nInd(1:N), ' -> 1D index:', ind
     print *, ' Check:(0?):', ind-i
     print *
  end do

end

