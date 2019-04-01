!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!c  Hij part of OSB: Only used for testing of Hij calculation               c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine getH0(CNT, IND, H0)
      integer, intent(IN)  :: cnt,  IND(cnt)
      double precision, intent(OUT) :: H0(cnt, cnt)

      integer :: I, J

      do I = 1, cnt
         do J = I + 1, cnt
           H0(I, J) = getHOSBHij(ind(I), ind(J))   
           H0(J, I) = H0(I, J)                 
         end do
         H0(I, I) =  getHOSBHIJ(ind(I), ind(I))
      end do

  end subroutine 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Calculate the specified HIJ = V^T|H|V  via HOSB        c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  double precision function getHOSBHij(row, col)      
      integer, intent(IN)  :: row, col

      integer :: level, bkIndex, rowIndex, colIndex
      double precision, allocatable :: Vi(:), Vj(:)
      integer :: TMP, nsize   

      if (row == col) then
          getHOSBHIJ = Eig0(row);  return
      end if
 
      call getLevelIndex0(row, col, level, bkIndex, rowIndex, colIndex)

      TMP = myHOSB%mStart(level)+((bkIndex-1)*sN(level)*sN(level)   &
            + (colIndex-1)*sN(level) + (rowIndex-1)) * myDim(level)

      if (level == 1) then
	  getHOSBHIJ = HOSB(tmp);  return
      end if

      nSize = myDim(level)

      allocate(VI(nSize),VJ(nSize))

      call getLevelV0(level-1, row, nSize, vi)
      call getLevelV0(level-1, col, nSize, vj)

      getHOSBHij=dot_product(Vi(1:nSize)*Vj(1:nSize),HOSB(TMP:tmp+nSize-1)) 

      deallocate(Vi)
      deallocate(Vj)

 end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

