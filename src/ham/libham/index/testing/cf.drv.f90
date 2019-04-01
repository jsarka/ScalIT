!
! Testing for CG and get3index 
!
PROGRAM TEST_CF_CFj
   implicit none

   integer :: jmax 
   INTEGER, allocatable :: jSize(:),cfInd(:,:),cfjInd(:,:)
   integer :: getCFSize, getCFPos, getCFjSize, getCFjPos
   integer :: ind1(4), ind2(2)
   integer :: inde0, inde1, js1,js2
   integer ::  I, OPT = 1

   DO WHILE (OPT /= 0)
       PRINT *, 'Input jmax CF/CFj '
       read  *, jmax
       allocate(jSize(jmax+1))
       print *, 'Input size from j=0 to j=',jmax,'. Total:',jmax+1
       read *, jSize(1:jmax+1)
  
       js1 = getCFSize(jmax,jSize); js2=getCFjSize(jmax,jSize)
       print *,'Size:CF=',js1,' CFj=',js2
       allocate(cfInd(4,js1), cfjInd(2,js2))
       call getCFIndex(jmax,jSize,js1,cfind)
       call getCFjIndex(jmax,jSize,js2,cfjind)

       print *
       print *, '# of Index for CF:', js1
       print *, 'Index, j1, j2,  j,  m1, K'
       DO I = 1, js1   
          print *, I, cfind(:,I)
       END DO

       print *
       print *, '# of index for CGFj:',js2
       print *, 'Index, j1, j2,  j,  m1, K'
       DO I = 1, js2 
          print *, I, cfjind(:,I)
       END DO

       deallocate(cfInd, cfjInd)

       print *
       print *, ' index<->(jm) for CF:'
       print *, 'Index   J   M   New Pos    Error'
       do inde0 = 1, js1   
          call getCFInd(jmax, jSize,inde0, ind1)
          inde1 = getCFPos(jmax,jSize, ind1)
          print *, inde0, ind1, inde1, ' DE:', (inde0-inde1)
       end do

       print *
       print *, ' index<->(jm) for CFj:'
       print *, 'Index   J   M   New Pos    Error'
       do inde0 = 1, js2
          call getCFjInd(jmax, jSize,inde0, ind2)
          inde1 = getCFjPos(jmax, jSize,ind2)
          print *, inde0, ind2, inde1, ' DE:', (inde0-inde1)
       end do
       deallocate(jSize)
       print *
       PRINT *, 'INPUT OPTION for loop: O FOR EXIT'
       READ *, OPT
    END DO
END
