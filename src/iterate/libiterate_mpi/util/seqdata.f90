!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                Auxillary subroutines for MOSBType                     c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine seqDataLenLong(gLen, nNodes, nID, pLen, pStart)
        implicit none
        include 'mpif.h'
        integer(KIND=MPI_OFFSET_KIND), intent(IN) :: gLen 
        integer, intent(IN)  :: nNodes, nID
        integer, intent(OUT) :: pLen
        integer(KIND=MPI_OFFSET_KIND), intent(OUT) :: pStart

        integer :: pNum

        pLen = gLen / nNodes
        pNum = gLen - pLen * nNodes

        if (nID<pNum) then
           pLen = pLen + 1
           pStart = nID * pLen + 1
        else
           pStart = pNum*(pLen+1) + (nID-pNum)*pLen + 1
        end if        

    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine seqDataPosLong(sF, gSize, gLen, gStart, gEnd)
        implicit none
        include 'mpif.h'
        integer, intent(IN) :: sF
        integer(kind=MPI_OFFSET_KIND), intent(IN) :: gSize(sF)
        integer(kind=MPI_OFFSET_KIND), intent(OUT):: gLen, gStart(sF), gEnd(sF)

        integer :: i
        gStart(1)=1; gEnd(1)=gSize(1)
        do i = 1, sF-1
           gStart(i+1) = gEnd(i)+1
           gEnd(i+1)   = gEnd(i)+gSize(i+1)
        end do
        gLen = gEnd(sF)

    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine seqDataPos(sF, pSize, pLen, pStart, pEnd)
        implicit none
        include 'mpif.h'
        integer, intent(IN) :: sF
        integer, intent(IN) :: pSize(sF)
        integer, intent(OUT):: pLen, pStart(sF), pEnd(sF)

        integer :: i
        pStart(1)=1; pEnd(1)=pSize(1)
        do i = 1, sF-1
           pStart(i+1) = pEnd(i)+1
           pEnd(i+1)   = pEnd(i)+pSize(i+1)
        end do
        pLen = pEnd(sF)

    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

