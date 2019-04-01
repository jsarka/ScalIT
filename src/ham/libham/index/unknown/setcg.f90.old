!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Subroutine to setup CG data: C(j1j2j;m(K-m)K)           c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine setCG(FNMax, jkNum, jkInd, mBase, jkmNum, mInd, cgdata)
   implicit none
   integer, intent(IN) :: FNMax, jkNum, jkmNum
   integer, intent(IN) :: jkInd(4, jkNum),mBase(jkNum),mInd(jkmNum)
   double precision, intent(OUT)  :: cgdata(jkmNum)

   double precision :: cg, LFN(FNMAX)
   integer :: i1, i2, j1, j2, j, m1, m2, m

   call lnFN(FNMAX, LFN)
   
   do i1 = 1, jkNum-1
      j1 = jkInd(1, i1);  j2 = jkInd(2, i1)
      j  = jkInd(3, i1);  m  = jkInd(4, i1)          
      do i2 = mBase(i1), mBase(i1+1)-1
         m1 = mInd(i2)
         m2 = m - m1
         cgdata(i2) = cg(j1, m1, j2, m2, j, m, FNMAX, LFN)
      end do    
   end do

   i1 = jkNum
   j1 = jkInd(1, i1);  j2 = jkInd(2, i1)
   j  = jkInd(3, i1);  m  = jkInd(4, i1)          
   do i2 = mBase(i1), jkmNum
      m1 = mInd(i2)
      m2 = m - m1
      cgdata(i2) = cg(j1, m1, j2, m2, j, m, FNMAX, LFN)
   end do

end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                         c
!c   Get the block data of CGdata and mindex. make sure    c
!c   that mcgNum = mSize(ind), otherwise there is messy    c
!c   Each block has the same (j1j2jK), but different m     c
!c                                                         c
!c jkNum: # of (j1,j2,j, K)                                c
!c jkInd: [4,jkNum], (j1,j2,j,K)                           c
!c mbase: [jknum], Base address for each (j1,j2,j,K) block c
!c jkmNum:# of (j1,j2,j; m,(K-m), K)                       c
!c mInd:  [jkmNum], the m index for (j1,j2,j; m,(K-m), K)  c
!c cgdata:[jkmNum], the CG value of C(j1j2j;m(K-m)K)       c
!c ind:   the index of the extracted (j1,j2,j,K) block     c
!c mcgNum: the number of m for the extracted (j1,j2,j,K)   c
!c jk:  [4], the extracted (j1,j2,j,K) indices             c
!c m0:  [mcgNum], the m indices for extracted (j1,j2,j,K)  c
!c cg0: [mcgNum], the corresponding CG(j1j2j;m(K-m)K)      c
!c                                                         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getJKM(jknum, jkInd, mBase, jkmNum, mInd, cgdata, ind, mcgNum, jk, m0, cg0)
   implicit none
   integer, intent(IN) :: jkNum, jkmNum, ind, mcgNum
   integer, intent(IN) :: jkInd(4,jkNum), mBase(jkNum), mInd(jkmNum)
   double precision,intent(IN) :: cgdata(jkmNum)

   integer, intent(OUT) :: jk(4), m0(mcgNum)
   double precision, intent(OUT) :: cg0(mcgNum)    

   jk(1:4) = jkInd(1:4, ind)
   m0(1:mcgNum) = mInd(mBase(ind):mBase(ind)+mcgNum-1)
   cg0(1:mcgNum)= cgdata(mBase(ind):mbase(ind)+mcgNum-1)

end

integer function stepSum(N, mSize, mSum)
   implicit none
   integer, intent(IN)  :: N, mSize(N)
   integer, intent(OUT) :: mSum(N)

   integer :: i

   mSum(1) = 1
   do i = 1, N-1
      mSum(i+1) = mSum(i)+mSize(i)
   end do
   stepSum = mSum(N)+mSize(N)
end
   


