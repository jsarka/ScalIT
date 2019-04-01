!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c              initialize MPI enviroment                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine preInit(ierr)     
     integer, intent(OUT) :: ierr

     logical :: mFlag
     double precision :: db
     double complex   :: cx

     call MPI_INITIALIZED(mFlag, ierr)
     if (.NOT. mFlag)  call MPI_INIT(ierr)
 
     call MPI_COMM_RANK(MPI_COMM_WORLD, id, ierr)
     call MPI_COMM_SIZE(MPI_COMM_WORLD, nNodes, ierr)

     call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, dbSize, ierr)
     call MPI_TYPE_SIZE(MPI_DOUBLE_COMPLEX, cxSize, ierr)

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine postInit(ierr)
     integer, intent(OUT) :: ierr

     integer :: i   

     call calGData(sF, sN, myconf)
     call calGDMNode(id, nNodes, myconf, myNode)
     call calMData(myconf, myNode, myData)
     call calViSeq(sF, sN, myVi)
     call calH0Seq(sNDVR, sF, sN, myH0)     
     call calHOSBGrid(sST, myconf, myData, myHOSB)
     call calVOSBGrid(myconf, myData, myVOSB)
     call calRESGrid(myconf, myData, myRES)
     call calDepGrid(sDep, myconf, myData, myDep)

     myNode%commID(sF) = MPI_COMM_WORLD
     do i = 1, sF-1
        if (myconf%gBlk(i)<nNodes) then
            call MPI_COMM_SPLIT(MPI_COMM_WORLD,myNode%grpID(i),   &
                            myNode%myID(i),myNode%commID(i),ierr)
        else
            myNode%commID(i) = MPI_COMM_SELF
        end if
     end do

!     do i = 1, sF-1
!        if (myNode%nodNum(i)>1) then
!            call MPI_COMM_SPLIT(MPI_COMM_WORLD,myNode%grpID(i),   &
!                            myNode%myID(i),myNode%commID(i),ierr)
!        else
!            myNode%commID(i) = MPI_COMM_SELF
!        end if
!     end do   
    
     blk(1:sF)  = myData%pBlk(1:sF)
     nin(1:sF)  = myData%pDim(1:sF)
     nout(1:sF) = sN(1:sF)
     pLen(1:sF) = blk(1:sF)*nin(1:sF)*nout(1:sF)

     pmax  = maxval(plen(1:sF))
     sNMax = maxval(sN(1:sF)) 

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function allocMOSB()
    integer :: info

    call deallocMOSB()

    allocMOSB=.false.
    
    allocate(VOSB(myVOSB%pLen), RES(myRES%pSize(sF)),        &
             EIG0(pmax), ResSeq(myRES%pSize(1)),stat=info)
    if (info/=0) return  

    select case (sJOB)
    case (JOB_RES1, JOB_RES2)
        allocate(AP(myRES%pSize(1)),stat=info)

    case (JOB_CRP,JOB_CRP1,JOB_CRP2)
        allocate(APP(myRES%pSize(1)), APR(myRES%pSize(1)),   &
                 AP(myRES%pSize(1)), stat=info)
    case default
        info=0
    end select

    if (info/=0) return

    if (sCX) then
       allocate(H0CX(myH0%mLen), stat=info)
       if (info/=0) return
       
       if (sST) then
         allocate(HOSBCX(myHOSB%pLen), stat=info)
         if (info/=0) return
       end if

       if (sNDVR) then
          allocate(OUTHCX(myHOSB%pSize(sF)),stat=info)
          if (info/=0) return
       end if

       if (totalDep) then
          allocate(DEPCX(myDEP%pLen), stat=info)
          if(info/=0) return
       end if
    else
       allocate(H0(myH0%mLen), stat=info)
       if (info/=0) return
       
       if (sST) then
         allocate(HOSB(myHOSB%pLen), stat=info)
         if (info/=0) return
       end if

       if (sNDVR) then
          allocate(OUTH(myHOSB%pSize(sF)),stat=info)
          if (info/=0) return
       end if

       if (totalDep) then
          allocate(DEP(myDEP%pLen),stat=info)
          if (info/=0) return
       end if
    end if

    allocMOSB = .true.
    
end function
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function allocHOSB(nSize)
   integer, intent(IN) :: nSize

   integer :: info

   allocHOSB = .TRUE.

   if (allocated(HOSB)) then

       if (size(HOSB) == nSize)  return
         
       deallocate(HOSB)

   end if

   allocate(HOSB(nSize), stat=info)

   allocHOSB= (info==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function allocHOSBCX(nSize)
   integer, intent(IN) :: nSize

   integer :: info

   allocHOSBCX = .TRUE. 
 
   if (allocated(HOSBCX)) then

       if (size(HOSBCX) == nSize) return

       deallocate(HOSBCX)

   end if
   
   allocate(HOSBCX(nSize), stat=info)

   allocHOSBCX= (info==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine destroyBJComm(ierr)
   integer, intent(OUT) :: ierr

   integer :: i

!   do i = 1, sF
!      if ((myNode%commID(i) /= MPI_COMM_WORLD) .AND.     &
!          (myNode%commID(i) /= MPI_COMM_SELF) )          &
!          call MPI_COMM_FREE(myNode%commID(i), ierr)
!   end do

    do i = 1, sF-1
       if (myconf%gBlk(i)<nNodes)    &
            call MPI_COMM_FREE(myNode%commID(i),ierr)
    end do

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine deallocMOSB()

    if (allocated(VOSB))   deallocate(VOSB)
    if (allocated(EIG0))   deallocate(EIG0)
    if (allocated(RES))    deallocate(RES)
    if (allocated(ResSeq)) deallocate(ResSeq)

    if (allocated(APP))    deallocate(APP)
    if (allocated(APR))    deallocate(APR)
    if (allocated(AP))     deallocate(AP)
   
    if (allocated(H0))     deallocate(H0)
    if (allocated(H0CX))   deallocate(H0CX)

    if (allocated(HOSB))   deallocate(HOSB)
    if (allocated(HOSBCX)) deallocate(HOSBCX)

    if (allocated(OUTH))   deallocate(OUTH)
    if (allocated(OUTHCX)) deallocate(OUTHCX)

    if (allocated(DEP))    deallocate(DEP)
    if (allocated(DEPCX))  deallocate(DEPCX)

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

