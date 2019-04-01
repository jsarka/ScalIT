!
! Testing for get4KSize and get4Kindex 
!
PROGRAM TEST_A2B2_A1K
   implicit none

   integer :: jTol,j1max, j2max, jmax, kmax
   logical :: parity, printIt
   INTEGER, allocatable :: jkmInd(:,:), kSize(:)
   integer :: get4kIndex,  OPT = 1, rSize
   integer :: i, j, k

   DO WHILE (OPT /= 0)
       PRINT *, 'Input Parity, JTOL, j1max, j2max, jmax, (T for even, F for odd parity)'
       read  *, parity, JTOL, j1max, j2max, jmax
       kmax = min(jmax, jtol)  
       allocate(kSize(kmax+1))

       call get4KSize(parity, JTOL, j1max, j2max, jmax, kmax, kSize)
       print *
       print *, 'kMax:', kmax
       print *, 'kSize:', kSize

       do i=0, kmax
          allocate(jkmInd(5,kSize(i+1)))
          rSize = get4kIndex(parity, jtol, j1max, j2max, jmax, i, kSize(i+1),jkmInd)
          print *
          print *, 'k:', i, '  kSize:', kSize(i+1), ' Real Size:',rSize
          print *, ' Print jkm Index:(T=Print, F=Not Print) '
          read(*, *) printIt
          if (printIt) then
             print *, 'jkm Index:'
             do j = 1, kSize(i+1)
                print *, jkmInd(1:5, j)
             end do
          end if
          deallocate(jkmInd)
       end do

       deallocate(kSize)
       PRINT *, 'INPUT OPTION for loop: O FOR EXIT'
       READ *, OPT
    END DO
END
