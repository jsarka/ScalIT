!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!c  Vi part of OSB:                                                       c
!c  getV(int col, double[n] v)  !V=V(f)*V(f-1)*...*V(1)                   c
!c  getLevelV(int level, int col, double[] v)                             c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      Get the blk and column indices for the specific column      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getLevelIndex(cnt, ind,  bkInd, colInd)
   integer (kind=8), intent(IN) :: cnt
   integer, intent(IN) :: Ind(cnt)
   integer, intent(OUT) :: bkInd(sF,cnt), colInd(sF,cnt)
   integer :: i, j 

   do i = 1, cnt
      do j = sF, 1, -1
         bkInd(j,i)  = (ind(i)-1)/myClen(j)
         colInd(j,i) = (ind(i)-bkInd(j,i)*myCLen(j)-1)/myDim(j) + 1          
      end do
   end do

   bkInd(1:sF,1:cnt) = bkInd(1:sF,1:cnt) + 1

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      Get the eigen-vector for the specific column        c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      VL(N,N) = V(level)*V(level-1)*...*V(2)*V(1)         c       
!c                       Vi = VL(:,col)                     c
!c   VL(F) = V, len=sN(level)*myDim(level)=myDim(level+1)   c
!c                 level = 0, 1, 2, ..., sF                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine getLevelVi(Level, bkInd, colInd, viLen, VI)
      integer, intent(IN) :: level, bkInd(sF),colInd(sF),viLen
      double precision, intent(OUT) :: VI(viLen)

      integer :: i, j, ind
      double precision :: tmp

      if (level==0) then
          Vi(1:viLen)=1.0D0;    return
      end if

      Vi(1:sN(1)) = 1.0D0
      
      do I = 1, level   
         do j = 2, sN(I)
            Vi((j-1)*myDim(i)+1 : j*myDim(i)) = Vi(1:myDim(i))
         end do

         ind = myVOSB%mStart(I)+(bkInd(I)-1)*sN(I)*sN(I)+    &
                  (colInd(i)-1)*sN(I)         

         do j = 1, sN(I)
              tmp = VOSB(ind+j-1)
              Vi((j-1)*myDim(i)+1:j*myDim(i)) =              &
                   tmp*Vi((j-1)*myDim(i)+1:j*myDim(i))
         end do
      end do

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      Get the eigen-vector for the specific column        c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   VA(N,N) = V(F)*V(F-1)*...*V(2)*V(1),  Vi = VA(:,col)   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine getVi(bkInd, colInd, VI)
      integer, intent(IN) :: bkInd(sF),colInd(sF)
      double precision, intent(OUT) :: VI(myLen)

      call getLevelVi(sF, bkInd, colInd, myLen, VI)

 end subroutine 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine getViVec(level,bkInd,colInd, allVi, oneVi)
      integer, intent(IN) :: level,bkInd,colInd
      double precision, intent(IN) :: allVi(sN(level),sN(level),myBlk(level))
      double precision, intent(OUT):: oneVi(sN(level)) 
      
      oneVi(1:sN(level)) = allVi(1:sN(level),colInd,bkInd)

 end subroutine 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
