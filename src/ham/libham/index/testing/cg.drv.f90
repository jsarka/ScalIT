!
! Testing for CG and get3index 
!
PROGRAM TEST_CG
   implicit none

   integer :: j1max, j2max, jmax,kmax
   INTEGER, allocatable :: jmkind(:,:)
   integer :: getCGSize, getCGPos, jmk(5)
   integer :: jsize, ind, ind1
   integer ::  I, OPT = 1

   DO WHILE (OPT /= 0)
       PRINT *, 'Input j1max,j2max,jmax,kmax for CG '
       read  *, j1max,j2max,jmax, kmax
  
       jsize = getCGSize(j1max,j2max,jmax,kmax)
       print *
       print *, '# of Index:', jsize

       allocate(jmkind(5,jSize))
       call getCGIndex(j1max,j2max,jmax,kmax,jSize,jmkind)
       print *, 'Index, j1, j2,  j,  m1, K'
       DO I = 1, jSize        
          print *, I, jmkind(:,I)
       END DO

       deallocate(jmkInd)

       print *, ' index<->(jm):'
       print *, 'Index   J   M   New Pos    Error'
       do ind = 1, jsize   
          call getCGInd(j1max,j2max,jmax,kmax, ind, jmk)
          ind1 = getCGPos(j1max,j2max,jmax,kmax, jmk)
          print *, ind, jmk, ind1, ' DE:', (ind-ind1)
       end do
       PRINT *, 'INPUT OPTION for loop: O FOR EXIT'
       READ *, OPT
    END DO
END
