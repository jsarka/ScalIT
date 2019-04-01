!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   The grouping of nodes for each block           c
!c     This only works when node > blk              c
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program test_NodeLen
  include 'mpif.h'
  integer(kind=MPI_OFFSET_KIND) :: gLen, gStart
  integer :: node, i, pLen, sLen, total, sstart, send

  print *, 'Input block Length and Node number: Length >> node'
  read (*,*) gLen, node
  
  print *, ' Call indexNodeLen:'
  do i = 0, node-1
     call indexNodeLen(gLen,node,i, pLen, gStart)
     write(*,100)  node, i, pLen, gStart
  end do

  print *, ' call MPI_NODE_LEN'
  total=gLen
  do i = 0, node-1 
    call getMPI1DIndex(total,node,i, sStart, sEnd, sLen)
     write(*,110)  node, i, sStart, sEnd, sLen
  end do

  100 format ('Total node:', I4, ' Node:',I4, ' : Plen:',I8, ' GStart:',I10)
  110 format ('Node #:', I4, ' Node:',I2, ' : Start:',I4, ' End:',I4, ' Len:',I4)

end

