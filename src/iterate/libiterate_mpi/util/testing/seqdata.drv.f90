!
! Program to test seqDataLen(), seqDataPos(), seqDataPosLong()
!
program test_seqdata
   implicit none
   include 'mpif.h'
   integer, parameter :: sMax=200, NMAX1=10000, NMAX2=900
   integer :: sF, id, nNode, i
   integer(kind=MPI_OFFSET_KIND) :: gLen,p0, gSize(sMAX), gStart(sMAX), gEnd(sMAX)
   integer :: pLen, pSize(sMAX), pStart(sMAX), pEnd(sMAX)
   double precision :: dM(sMAX)
   integer :: opt = 1

   print *, '==================================='
   print *, '     Testing seqData subroutines'
   print *, '===================================='
   
   do while (opt/=0)
      print *, 'input sF and nNode'
      read(*,*) sF, nNode
   
      call random_number(dM)
      gSize(1:sF) = NMAX1*dM(1:sF)
      pSize(1:sF) = NMAX2*dM(1:sF)
      
      do i = 1, sF
         if (gSize(i)==0) gSize(i)=NMAX1/2
         if (pSize(i)==0) pSize(i)=NMAX2/2
      end do

      gLen = gSize(1)

      print *, ' Distribute the data:'
      print *, 'gLen=', gLen, ' Nodes=', nNode

      do id = 0, nNode-1
         call seqDataLenLong(gLen, nNode, id, pLen, p0)
         print *, 'id=', id, ' pLen=', pLen, ' Pstart=', p0
      end do

      print *
      print *, ' Global data index:'
      call seqDataPosLong(sF, gSize, gLen, gStart, gEnd)

      print *, ' Total Length:', gLen
      print *, '       i        gSize     gStart      gEnd'
      do i = 1, sF
         print *, i, gSize(i), gStart(i), gEnd(i)
      end do

      print *
      print *, ' Local data index:'
      call seqDataPos(sF, pSize, pLen, pStart, pEnd)

      print *, ' Total Length:', pLen
      print *, '        i        pSize     pStart      pEnd'
      do i = 1, sF
         print *, i, pSize(i), pStart(i), pEnd(i)
      end do
   
      print *, 'Input options for continue: 0=exit, other=continue'
      read(*,*) opt

   end do

   print *
   print *, ' -----------   Finish Testing -------------- '

end program
