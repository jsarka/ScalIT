!
! Testing subroutines in index_nlm   
!
PROGRAM TEST_index_nlm
   implicit none

   integer :: nmax
   INTEGER, allocatable :: nlmInd(:,:)
   integer :: getNLMSize, getNLMPos
   integer :: nsize, ind, ind1, N, L, m, I, OPT = 1

   DO WHILE (OPT /= 0)
       PRINT *, 'Input Nmax for NLM '
       read  *, nmax
  
       nsize = getNLMSize(nmax)
       print *
       print *, '# of Index:', nsize

       allocate(nlmInd(3,nSize))
       call getNLMIndices(nSize, nlmInd)
       print *, 'Index,      nIndex,       LIndex       mIndex'
       DO I = 1, nSize        
          print *, I, nlmInd(1:3,I)
       END DO

       deallocate(nlmInd)

       print *, ' index<->(NLM):'
       print *, 'Index           N            L          M         New Pos    Error'
       do ind = 1, nsize   
          call getNLMIndex(ind, n, L, m)
          ind1 = getNLMPos(N, L, m)
          print *, ind, N, L, M, ind1, 'DE:', (ind-ind1)
       end do
       PRINT *, 'INPUT OPTION for loop: O FOR EXIT'
       READ *, OPT
    END DO
END
