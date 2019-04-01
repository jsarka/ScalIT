!cccccccccccccccccccccccccccccccccccccc
!c   test program for get***Index     c
!cccccccccccccccccccccccccccccccccccccc
program test_getvindex
   implicit none
   include 'mpif.h'
   integer, parameter :: NMAX=20
   integer :: sF, sN(NMAX)
   integer(kind=MPI_OFFSET_KIND) :: blkInd1(NMAX),blkInd2(NMAX)
   integer :: sNInd1(NMAX),snInd2(NMAX)   
   integer(kind=MPI_OFFSET_KIND) :: nlen, pos, pos1,colInd(NMAX) !, getViPos

   integer :: i, j

   read(*,*) sF
   read(*,*) sN(1:sF)

   print *
   print *, ' Layer:',sF
   print *, ' Layer Configuration:', sN(1:sF)
   print *, ' If there is no output, everything is OK!'

   nlen=1
   do i = 1, sF
       nLen = nLen*sN(i)
   end do

   do i = 1, nlen
!      i = nlen
!      print *
!      print *, ' Position=',i
      pos = i
      call getViColInd(sF,sN,pos,blkInd1,snInd1)
      call getViPosInd(sF,sN,pos,blkInd2,snInd2,colInd)

      do j = 1, sF
         if ((blkInd1(j)/=blkInd2(j)) .OR. (snInd1(j)/=snInd2(j)))   &
            print *, 'Error in getVi(Col/Pos)Ind at level=',j,       &
                      blkInd1(j),blkInd2(j),snInd1(j),snInd2(j)
        ! else
        !    print *, ' GetVi(Col/Pos)Ind, getVi(Col/Pos)Index are OK!'
        ! end if

         call getViPos(sF,sN,j,blkInd2(j),sNInd2(j),colInd(j),pos1)
         if (pos1 /= pos)     &
             print *, 'Error in getViPos:  Level=',j,blkInd2(j),snInd2(j),colInd(j),pos1  
        ! else
        !     print *, ' GetViPos, getViPosition are OK!'
        ! end if 
      end do

   end do

end program
