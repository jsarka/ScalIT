!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   The grouping of nodes for each block           c
!c     This only works when node > blk              c
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program test_NodeLen

  include 'mpif.h'
  integer(kind=MPI_OFFSET_KIND) :: gLen, gStart
  integer :: node, i, pLen

  print *, 'Input block Length and Node number: Length >> node'
  read (*,*) gLen, node
  
  do i = 0, node-1
     call indexNodeLen(gLen,node,i, pLen, gStart)
     write(*,100)  node, i, pLen, gStart
  end do

  100 format ('Total node:', I4, ' Node:',I4, ' : Plen:',I8, ' GStart:',I10)

end

