!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   The grouping of nodes for each block           c
!c     This only works when node > blk              c
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program test_node

  integer :: blk, node, i, grp, id, gpSize

  print *, 'Input block Number and Node number: block < node'
  read (*,*) blk, node
  
  do i = 0, node-1
     call indexNode(blk,node,i, id, grp, gpSize)
     write(*,100)  i, grp, id, gpSize
  end do

  100 format ('node:', I4, ' GroupID:',I4, ' : ID in subgroup:',I4, ' Group Size:',I4)

end

