!
! Testing for geti4Size and get4index 
!
PROGRAM TEST_ABCD
   implicit none
   integer, parameter :: NUM = 100

   integer :: jTol,j1max, j2max, jmax
   logical :: parity
   INTEGER, allocatable :: jkInd(:,:), mind(:), msizeind(:)
   integer :: jsize,msize, jsize1, j2size(2), I, OPT = 1

   DO WHILE (OPT /= 0)
       PRINT *, 'Input Parity, JTOL, j1max, j2max, jmax, (T for even, F for odd parity)'
       read  *, parity, JTOL, j1max, j2max, jmax
  
       call get4Size(parity, JTOL, j1max, j2max, jmax, jsize, msize)
       print *
       print *, '# of (j1j2jK):', jsize, ' # of CG:',mSize

       allocate(jKInd(4,jSize),msizeind(jsize),mind(msize))
       call get4Index(parity, jtol, j1max, j2max, jmax, jSize, jkInd, msizeind, msize, mInd)
       print *, 'Index, j1Index, j2index, jindex, kIndex'
       jsize1 = min(jsize, NUM)
       DO I = 1, jSize1
          print *, I, jKInd(1:4,I), msizeind(I)
       END DO

       print *
       print *, 'mIndex:'
       print *,  mind(1:jsize1)

       deallocate(jKInd, msizeind, mind)
       PRINT *, 'INPUT OPTION for loop: O FOR EXIT'
       READ *, OPT
    END DO
END
