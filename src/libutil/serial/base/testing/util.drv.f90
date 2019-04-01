!
! Testing subroutine  for Util_cx
!

program util_drv
   implicit none
   integer, parameter :: N = 10
   double precision, dimension(N) :: vec, oldVec, newVec   
   double precision :: minval, maxVal,  val, E0
   double precision :: getMin, getMax, getMaxUnder, getMinAbove
   double precision :: getMin0, getMax0, getMaxUnder0, getMinAbove0
   integer, dimension(N) :: newIndex
   integer :: minInd, maxInd, ind, I, M
   integer :: ntype

   call random_number(vec)
   oldVec(1:N) = vec(1:N)
   print *, 'Original Vector:'
   call printVec(N, vec)
    
   print *
   print *, '***************************'
   print *, 'Min/Max Testing:'
   print *, '***************************'
   call getMinMax(N, vec, minVal, maxVal)
   print *, 'getMinMax: Min=', minVal, '  Max=', maxVal
   call getMinMax0(N, vec, minVal, minInd, maxVal, maxInd)
   print *, 'getMinMax0:Min=', minVal, '  Max=', maxVal
   print *, 'MinInd=', minInd, '  MaxInd=', maxInd
   minVal = getMin(N, vec)
   maxVal = getMax(N, vec)
   print *, 'getMin/getMax:Min=', minVal, '  Max=', maxVal
   minVal = getMin0(N, vec, minInd)
   maxVal = getMax0(N, vec, maxInd)
   print *, 'getMin0/Max0: Min=', minVal, '  Max=', maxVal
   print *, 'MinInd=', minInd, '  MaxInd=', maxInd

   print *, '***************************'
   print *, '     Real Ascend/Descend   '
   print *, '***************************'  
   print *, 'Real Ascend Vector:'
   call Reorder('A', N, vec)
   call printVec(N, vec)
   
   vec(1:N) = oldVec(1:N)  
   print *, 'Real Descend Vector:'
   call Reorder('D', N, vec)
   call printVec(N, vec)


   print *, '***************************'
   print *, '    Test getNextMin(Max)'
   print *, '***************************'
   vec(1:N) = oldVec(1:N)
   print *, 'Descending' 
   maxVal = getMax0(N, vec, maxInd)
   print *, maxVal, oldVec(maxInd)
   DO I = 1, N-2
      val = getMaxUnder0(N, Vec, maxVal, ind)
      print *, val, oldVec(ind)
      maxVal = val
   END DO
   minVal = getMin0(N, vec, minInd)
   print *, minVal, oldVec(minInd)

   vec(1:N) = oldVec(1:N)
   print *, 'Ascending' 
   minVal = getMin0(N, vec, minInd)
   print *, minVal,  oldVec(minInd)
   DO I = 1, N-2
      val = getMinAbove0(N, Vec, minVal, ind )
      print *, val,  oldVec(ind)
      minVal = val
   END DO
   maxVal = getMax0(N, vec, maxInd)
   print *, maxVal,  oldVec(maxInd)

!******************************
   vec(1:N) = oldVec(1:N)
   minVal = getMin(N, vec)
   E0 = minVal
   M=N
   print *
   print *, '*************  Test GetWindow ***********************'
   print *, 'Energy Window: E0=', E0 , 'Window Number:',M
   CALL getWindow(E0, N, vec, M, newVec)
   call  printVec(M, newVec)  


   vec(1:N) = oldVec(1:N)
   minVal = getMin(N, vec)
   E0 = minVal
   M=N/2
   print *
   print *, 'Energy Window: E0=', E0 , 'Window Number:',M
   CALL getWindow(E0, N, vec, M, newVec)
   call  printVec(M, newVec)  

   vec(1:N) = oldVec(1:N)
   maxVal = getMax(N, vec)
   E0 = maxVal
   M=N
   print *
   print *, 'Energy Window: E0=', E0 , 'Window Number:',M
   CALL getWindow(E0, N, vec, M, newVec)
   call  printVec(M, newVec)  

   vec(1:N) = oldVec(1:N)
   maxVal = getMax(N, vec)
   E0 = maxVal
   M=N/2
   print *
   print *, 'Energy Window: E0=', E0 , 'Window Number:',M
   CALL getWindow(E0, N, vec,M, newVec)
   call  printVec(M, newVec)  


   vec(1:N) = oldVec(1:N)
!   minVal = getMin(N, vec, minInd)
   E0 = (minVal+maxVal)/2.0D0 
   M=N
   print *
   print *, 'Energy Window: E0=', E0 , 'Window Number:',M
   CALL getWindow(E0, N, vec, M, newVec)
   call  printVec(M, newVec)  

   vec(1:N) = oldVec(1:N)
!   minVal = getMin(N, vec, minInd)
   E0 = (minVal+maxVal)/2.0D0 
   M=N/2
   print *
   print *, 'Energy Window: E0=', E0 , 'Window Number:',M
   CALL getWindow(E0, N, vec, M, newVec)
   call  printVec(M, newVec)  

end

