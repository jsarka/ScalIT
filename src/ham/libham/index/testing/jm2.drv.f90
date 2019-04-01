!
! Testing subroutines in index_jm2   
!
PROGRAM TEST_index_jm2
   implicit none

   integer :: jmax
   INTEGER, allocatable :: jmInd(:,:)
   integer :: getJm2Size, getJm2Pos
   integer :: nsize, ind, ind1, j, m, I, OPT = 1

   DO WHILE (OPT /= 0)
       PRINT *, 'Input jmax for jm1 '
       read  *, jmax
  
       nsize = getJm2Size(jmax)
       print *
       print *, '# of Index:', nsize

       allocate(jmInd(2,nSize))
       call getJm2Indices(nSize, jmInd)
       print *, 'Index,      jIndex       mIndex'
       DO I = 1, nSize        
          print *, I, jmInd(1:2,I)
       END DO

       deallocate(jmInd)

       print *, ' index<->(jM):'
       print *, 'Index           j          M         New Pos    Error'
       do ind = 1, nsize   
          call getJm2Index(ind, j, m)
          ind1 = getJm2Pos(j, m)
          print *, ind, j, M, ind1, 'DE:', (ind-ind1)
       end do
       PRINT *, 'INPUT OPTION for loop: O FOR EXIT'
       READ *, OPT
    END DO
END
