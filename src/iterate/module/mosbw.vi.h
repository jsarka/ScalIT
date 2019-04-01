!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!c  Vi part of OSB:                                                          c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Calculate Vi at specified level. Vi is distributed               c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getLevelVi(level, ind, Vi)
   integer, intent(IN) :: level, ind
   double precision, intent(OUT) :: Vi(myVi%mEnd(level))

   integer :: req(sF), reqNum, status(MPI_STATUS_SIZE), ierr
   integer :: i, Btag
   double precision :: tmpV(sNMax)

   reqNum = 0
   if (rootInd(level,ind)==id) then    ! receive data
      do i = 1, level
  
         if (rootInd(i,ind)==id) then   ! local copy
  	     call copy3LVec(sN(i),blk(i),VOSB(myVOSB%pStart(i)),blkInd(i,ind),&
                           sNInd(i,ind),Vi(myVi%mStart(i)))
         else       ! receive data
             reqNUM = reqNum + 1
  	     BTag = B0Tag*i +ind ; 
             call MPI_IRecv(Vi(myVi%mStart(i)),sN(i),MPI_DOUBLE_PRECISION,    &
                    rootInd(i,ind), BTag, MPI_COMM_WORLD, req(reqNum), ierr)
         end if
      end do

      do i = 1, reqNum
         call MPI_WAIT(req(i),status, ierr)
      end do

   else
      do i = 1, level-1
         if (rootInd(i,ind)==id) then
	    call copy3LVec(sN(i),blk(i),VOSB(myVOSB%pStart(i)),blkInd(i,ind), &
                           sNInd(i,ind),tmpV)
            BTag = B0Tag * i + ind
            call MPI_Send(tmpV,sN(i),MPI_DOUBLE_PRECISION,rootInd(level,ind), &
                          BTag, MPI_COMM_WORLD, ierr)
         end if
      end do

   end if

   if ((myNode%grpID(level)==grpInd(level,ind)).AND.(nodeNum(level,ind)>1)) then
       call MPI_BCAST(Vi,myVi%mEnd(level),MPI_DOUBLE_PRECISION,              &
                rootInd(level,ind), myNode%commID(level), ierr)
   end if

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getFullVi(ind, Vi)
   integer :: ind
   double precision, intent(OUT) :: Vi(myVi%mLen)

   integer :: req(FMAX), reqNum, status(MPI_STATUS_SIZE), ierr
   integer :: i, BTag
   double precision ::  tmpV(sNMax)
 
   reqNum = 0
   if (id==rootID) then          
      do i = 1, sF
         if (rootInd(i,ind)==id) then
            call copy3LVec(sN(i),blk(i),Vosb(myVosb%pStart(i)),blkInd(i,ind),&
 	           sNInd(i,ind),Vi(myVi%mStart(i):myVi%mEnd(i)))
         else
            reqNum = reqNum + 1
            BTag = i*B0Tag + ind
            call MPI_IRecv(Vi(myVi%mStart(i)),sN(i),MPI_DOUBLE_PRECISION,    &
                   rootInd(i,ind),BTag,MPI_COMM_WORLD,req(reqNum),ierr)
         end if
      end do
     
      do i = 1, reqNum
         call MPI_WAIT(req(i),status,ierr)
      end do
   else
      do i = 1, sF
         if (rootInd(i,ind)==id) then
            BTag = i*B0Tag + ind
            call copy3LVec(sN(i),blk(i),Vosb(myVosb%pStart(i)),blkInd(i,ind),&
 	           sNInd(i,ind),tmpV)
            call MPI_Send(tmpV,sN(i), MPI_DOUBLE_PRECISION, rootID, BTag,    &
                   MPI_COMM_WORLD, ierr)
         end if
      end do
   end if

   call MPI_BCAST(Vi,myVi%mLen,MPI_DOUBLE_PRECISION,rootID,MPI_COMM_WORLD,ierr)

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


