!
! Testing for get3Size and get3index 
!
PROGRAM TEST_ABC
   implicit none

   integer :: jTol,jmax
   logical :: parity, showit
   INTEGER, allocatable :: jInd(:), kInd(:)
   integer :: get3Size, jsize,I, OPT = 1

   DO WHILE (OPT /= 0)
       PRINT *, 'Input JTOL, jmax, parity(T for even, F for odd)'
       read  *, JTOL, jmax, parity
  
       jsize = get3Size(parity, JTOL, jmax)
       print *
       print *, '# of Index:', jsize

       print *, " Show (jk) index: T=Show it, F=Don't show it"
       read *, showit
       if (showit) then
          allocate(jInd(jSize),kInd(jSize))
          call get3Index(parity, jtol, jmax, jSize, jInd, kInd)
          print *, 'Index, jIndex, kIndex'
          DO I = 1, jSize        
             print *, I, jInd(I), kInd(I)
          END DO

          deallocate(jInd, kInd)
       end if

       PRINT *, 'INPUT OPTION for loop: O FOR EXIT'
       READ *, OPT
    END DO
END
