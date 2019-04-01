!
! Testing subroutine  for Util_cx
!

program util_cx_drv
   integer, parameter :: N = 10
   double complex, dimension(N) :: vec, oldVec   
   integer :: ntype

   call randVec_cx(N, vec)
   oldVec(1:N) = vec(1:N)
   print *, 'Original Vector:'
   call printVec_cx(N, vec)
  
   nType = 1
   print *
   print *, 'Real Ascend Vector:'
   call AscReorder_cx(nType, N, vec)
   call printVec_cx(N, vec)

   vec(1:N) = oldVec(1:N)
   print *
   print *, 'Real Descend Vector:'
   call DesReorder_cx(nType, N, vec)
   call printVec_cx(N, vec)

   nType = 3
   vec(1:N) = oldVec(1:N)
   print *
   print *, 'Imag Ascend Vector:'
   call AscReorder_cx(nType, N, vec)
   call printVec_cx(N, vec)

   vec(1:N) = oldVec(1:N)
   print *
   print *, 'Imag Descend Vector:'
   call DesReorder_cx(nType, N, vec)
   call printVec_cx(N, vec)

end

