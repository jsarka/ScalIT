!
! Testing for get3Size and get3index 
!
PROGRAM TEST_ABC
   implicit none

   integer :: jmax
   INTEGER, allocatable :: jkInd(:,:)
   integer :: getPjmSize, getPjmPos
   integer :: jsize, ind, ind1, j, m, I, OPT = 1

   DO WHILE (OPT /= 0)
       PRINT *, 'Input Jmax for Pjm '
       read  *, Jmax
  
       jsize = getPjmSize(jmax)
       print *
       print *, '# of Index:', jsize

       allocate(jKInd(2,jSize))
       call getPjmIndex(jmax, jkInd)
       print *, 'Index, jIndex, mIndex'
       DO I = 1, jSize        
          print *, I, jKInd(1:2,I)
       END DO

       deallocate(jKInd)

       print *, ' index<->(jm):'
       print *, 'Index   J   M   New Pos    Error'
       do ind = 1, jsize   
          call getPjmInd(jmax, ind, j, m)
          ind1 = getPjmPos(jmax, j, m)
          print *, ind, j, m, ind1, ' DE:', (ind-ind1)
       end do
       PRINT *, 'INPUT OPTION for loop: O FOR EXIT'
       READ *, OPT
    END DO
END
