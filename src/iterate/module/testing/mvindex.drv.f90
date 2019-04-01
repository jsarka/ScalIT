!
! Testing mindex.f90
!
program test_mindex
   use mosbtype
   implicit none

   include 'mpif.h'
   integer :: nNodes, sF, sN(FMAX)
   integer  :: i, j
   integer(kind=MPI_OFFSET_KIND) :: pos, pos1
   TYPE(GDataInfo) :: myData
   integer :: blkInd1(FMAX),snInd1(FMAX),colInd(FMAX)
   integer :: grpInd1(FMAX),rootInd1(FMAX), nodNum1(FMAX), idInd1(FMAX)

   integer :: blkInd2(FMAX),snInd2(FMAX),colInd2(FMAX)
   integer :: grpInd2(FMAX),rootInd2(FMAX), nodNum2(FMAX), idInd2(FMAX)

   read(*,*) nNodes, sF
   read(*,*) sN(1:sF)

   call calGData(sF,sN,myData)
   call printGData(myData)

   do j = 1, sF
      print *, ' Layer=',j
      print *, '-----------------------------------------------------------------------'
      print *, '  Pos           blkInd     sNInd    grpInd    rootInd    nNum '
      print *, '-----------------------------------------------------------------------'
      do i = 1, myData%gN
         pos = i
         call mgetViColInd(nNodes,sF,sN,pos,grpInd1,rootInd1,nodNum1,blkInd1,snInd1)     
         call mgetViPosInd(nNodes,sF,sN,pos,grpInd2,rootInd2,nodNum2,idInd2,blkInd2,snInd2,colInd2)
         call mgetViPos(j,nNodes,sF,sN,grpInd2(j),rootInd2(j),nodNum2(j),idInd2(j),blkInd2(j),   &
                        snInd2(j),colInd2(j),pos1)
         if ( (grpInd1(j)/=grpInd2(j)) .or. (rootInd1(j)/=rootInd2(j)) .or.   &
              (nodNum1(j)/=nodNum2(j)) .or. (blkInd1(j)/=blkInd2(j))   .or.   &
              (snInd1(j)/=snInd2(j)) )         &
            print *, 'D(grpInd)=',grpInd1(j)-grpInd2(j),' D(rootInd)=',       &
                   (rootInd1(j)-rootInd2(j)), ' D(nodNum)=',     &
                   (nodNum1(j)-nodNum2(j)), ' D(blkInd)=', (blkInd1(j)-blkInd2(j)), & 
                   ' D(snInd)=', (snInd1(j)-snInd2(j)) 
         if (pos1/=pos) print *, ' Diff of Position:',(pos1-pos), ' Pos=',pos, 'Pos1=',pos1

       ! write(*,10) pos,blkInd1(j),sNInd1(j),grpInd1(j),rootInd1(j),nodNum1(j), pos1-pos
!       if (j==2)  write(*,10) pos,blkInd2(j),sNInd2(j),grpInd2(j),rootInd2(j),colInd2(j),nodNum2(j)
!         write(*,10) pos,sN(j),blkInd1(j)-blkInd2(j),sNInd1(j)-sNInd2(j),   &
!                     rootInd1(j)-rootInd2(j),nodNum1(j)-nodNum2(j), pos1
      end do
   end do

   10 format (I15, 5I10, I6)

end
