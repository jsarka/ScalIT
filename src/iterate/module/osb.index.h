!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!c  Grid part of OSB_BASE:                                             c
!c                                                                     c
!c  get1D(): get 1D index from a MD index                              c
!c  getMD(): get the MD index fromm 1D index                           c
!c  getBLKD(): used to calculate V                                     c
!c  getRCIND():get block, rol, and column at each level                c 
!c  getLevelIndex(): get the level index, row, col, blk indices from   c
!c                   (row, col). Used to calculate Hij using HOSB      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Subroutines to convert index between 1D and mutidimension.   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccc
!c   Get 1D index from multidimension index    c
!ccccccccccccccccccccccccccccccccccccccccccccccc 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function get1D(MDIND)
     integer, intent(IN) :: MDIND(sF)

      integer :: I

      get1D = MDIND(1)
      do I = 2, sF
         get1D = GET1D + (MDIND(I)-1) * myDim(I)
      end do

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccc
!c       Get MD index from 1D index            c
!ccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine getMD(IND, MDIND)
      integer, intent(IN) :: IND
      integer, intent(OUT) :: MDIND(sF)

      integer :: TMP, I

      TMP = IND
      do I = sF, 2, -1
         MDIND(I) = (TMP-1) / myDim(I) + 1
         TMP = TMP - (MDIND(I)-1) * myDim(I)
      end do
      MDIND(1) = TMP

 end subroutine 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
!c   Get block index and column index in that block    c
!c   for the given column, used to find VOSB position  c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine getBlkd(col, BKIND,INBKD)    
      integer, intent(IN)               :: col
      integer, dimension(sF),intent(OUT) :: BKIND, INBKD

      integer :: I

      do I = 1, sF
         BKIND(i) = (col-1) / (myDim(i)*sN(i)) + 1
         INBKD(i) = (col-1) - ( BKIND(i) - 1 ) * (myDim(i)*sN(i))
         INBKD(i) = INBKD(i) / myDim(i) + 1
      end do

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
!c   Get level, block index and row and column index in the block     c
!c   for given point (row,col), used to cal HO(I, J) by using HOSB    c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine getLevelIndex0 (row, col, level, bkInd, rowInd, colInd)
    
      integer, intent(IN)               :: row, col
              ! which block the (row, col) point is in
              ! if bkind=0, not belong to any block
      integer, intent(OUT)  :: LEVEL   
              ! I > LEVEL, ROWIND=COLIND, I=LEVEL, ROWIND !=COLIND   
      integer, intent(OUT)  :: BKIND  
              ! block number at level, used for HOSB
      integer, intent(OUT)  :: ROWIND, COLIND  
              ! ROW AND COL INDICES IN THE BLOCK at level-1

      integer  :: I

      do I = sF, 1, -1
        rowInd  = (row-1) / myDim(I) + 1
        colInd  = (col-1) / myDim(I) + 1 
        bkInd   = (row-1) / (myDim(I)*sN(I)) + 1

        if ((rowInd /= colInd) .or. (I == 1) )  then
               level  = i

              rowInd = row - (bkInd-1) * sN(I) * myDim(I)
              rowInd =  (rowInd-1)/myDim(I) + 1

              colInd = col - (bkInd-1) * sN(I) * myDim(I) 
              colInd = (colInd-1)/myDim(I) + 1

              exit
        end if
      end do

 end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
!c   Get block index and row and column index in that block    c
!c   for the given point (row,col), used to cal HO(I, J)       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine getRCInd(row, col, BKIND, rowInd, colInd)
      integer, intent(IN)               :: row, col
                ! which block the (row, col) point is in, 
                ! if bkind=0, not belong to any block
      integer, dimension(sF),intent(OUT) :: BKIND
                ! row, and col index within the block, 
                ! rowInd(level) <= myDim(level)
      integer, dimension(sF),intent(OUT) :: rowInd, colInd  

      integer :: level

      BKIND(sF)  = 1                  
      rowInd(sF) = (row-1) / myDim(sF) + 1
      colInd(sF) = (col-1) / myDim(sF) + 1

      do level = sF-1, 1, -1         
         if (rowInd(level+1)  .ne. colInd(level+1))  then
            bkInd(level) = 0    ! no contribution at this level
            rowInd(level) = rowInd(level+1)
            colInd(level) = colInd(level+1)
         else
            bkInd(level) = (row-1) / (sN(level)*myDim(level)) + 1
            rowInd(level) = rowInd(level) - (bkInd(level)-1)     &
                            * sN(level) * myDim(level)
            colInd(level) = colInd(level) - (bkInd(level)-1)     &
                            * sN(level) * myDim(level)
         end if
      end do

      end subroutine GETRCIND
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
!c   Get level, block index and row and column index in the block     c
!c   for given point (row,col), used to cal HO(I, J) by using HOSB    c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getAllLevelIndex(cnt,rcIndex,levelInd,bkIndex,rowIndex,colIndex)
      integer (kind=8), intent(IN)           :: cnt
      integer,intent(IN), dimension(CNT)     :: rcIndex 
      integer,intent(OUT),dimension(CNT,CNT) :: levelInd, bkIndex,   &
                                                rowIndex, colIndex
      
      integer  :: I, JJ, KK
      integer  :: row, col, level, bkInd, rowInd, colInd

      do JJ = 1, CNT
	 do KK = 1, CNT
	    row = rcIndex(JJ)
            col = rcIndex(KK)
            do I = sF, 1, -1
               rowInd  = (row-1) / myDim(I) + 1
               colInd  = (col-1) / myDim(I) + 1 
               bkInd   = (row-1) / (myDim(I)*sN(I)) + 1

               if ((rowInd /= colInd) .or. (I == 1) )  then
                   level  = i

                   rowInd = row - (bkInd-1) * sN(I) * myDim(I)
                   rowInd =  (rowInd-1)/myDim(I) + 1

                   colInd = col - (bkInd-1) * sN(I) * myDim(I) 
                   colInd = (colInd-1)/myDim(I) + 1

                   exit
               end if
            end do
	    levelInd(JJ, KK) = level
            bkIndex (JJ, KK) = bkInd
            rowIndex(JJ, KK) = rowInd
	    colIndex(JJ, KK) = rowInd
         end do
      end do
end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

