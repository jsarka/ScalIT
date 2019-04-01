!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!c  Vi part of OSB:                                         c
!c  Subroutines:                                            c
!c  getV0(int col, double[n] v)  !V=V(f)*V(f-1)*...*V(1)    c
!c  getLevelV0(int level, int col, double[] v)              c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      Get the eigen-vector for the specific column        c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      VL(N,N) = V(level)*V(level-1)*...*V(2)*V(1)         c       
!c                       Vi = VL(:,col)                     c
!c      VL(F) = V , N=nSize=myDim(level)*sN(level)          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getLevelV0(Level, col, nSize, VI)
      integer, intent(IN) :: level
      integer, intent(IN) :: col, nSize   
      double precision, intent(INOUT) :: Vi(nSize)

      integer, dimension(sF) :: BIBLK, VINBK
      integer  :: i, j, ind
      double precision :: tmp

      call getBlKD(col, BIBLK, VINBK)

      if (level==0) then
         Vi(1:nSize)=1.0D0; return
      end if

      Vi(1:sN(1)) = 1.0D0
      
      do I = 1, level         
         do j = 2, sN(I)
            Vi((j-1)*myDim(i)+1:j*myDim(i))=Vi(1:myDim(i))
         end do

         ind = myVOSB%mStart(I)+(BIBLK(I)-1)*sN(I)*sN(I)+    &
                  (VINBK(i)-1)*sN(I)
         do j = 1, sN(I)
              tmp = VOSB(ind+j-1)
              VI((j-1)*myDim(i)+1:j*myDim(i)) =              &
                       tmp*VI((j-1)*myDim(i)+1:j*myDim(i))
         end do

      end do

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      Get the eigen-vector for the specific column        c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   VA(N,N) = V(F)*V(F-1)*...*V(2)*V(1),  Vi = VA(:,col)   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine getV0(col, VI)
      integer, intent(IN) :: col
      double precision, intent(INOUT) :: VI(myLen)

      call getLevelV0(sF, col, mylen, VI)

  end subroutine 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     
