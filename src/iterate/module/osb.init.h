!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calOSB()
   integer :: i

   myDim(1) = 1 
   do i = 1,sF-1
      myDim(i+1)=myDim(i)*sN(i)
   end do
 
   myBlk(sF) = 1
   do i = sF, 2, -1
      myBlk(i-1) = myBlk(i)*sN(i)
   end do

   myCLen(1:sF) = myDim(1:sF)*sN(1:sF)
   myLen = myBlk(1)*sN(1)

   outLen=0
   if (sNDVR)   outLen = sN(sF)*myLen

   myVi%mSize(1:sF) = sN(1:sF)
   call calDataInfo(sF, myVi)

   myH0%mSize(1:sF)   = sN(1:sF)**2  
   if (sNDVR) myH0%mSize(sF)=0
   call calDataInfo(sF, myH0)
   call printDataInfo(myH0)

   myHOSB%mSize(1:sF)    = sN(1:sF)*myLen
   call calDataInfo(sF, myHOSB)
   if (.NOT.sST) then
      myHOSB%mStart(1:sF) = 1
      myHOSB%mEnd(1:sF) = myHOSB%mSize(1:sF)
      myHOSB%mLen = myHOSB%mMaxSize
   end if

   myVOSB%mSize(1:sF) = myBlk(1:sF)*(sN(1:sF)**2)
   call calDataInfo(sF, myVOSB)

   myDep%mSize(1:sF)=0
   totalDep=.false.
   do i = 1, sF
     if (sDep(i)) then
       totalDep=.true.
       myDep%mSize(i)=myDim(i)
     end if
   end do
   call calDataInfo(sF, myDEP)

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function allocHOSB(nSize)
   integer, intent(IN) :: nSize

   integer :: info

   allocHOSB = .TRUE.

   if (allocated(HOSB)) then
       if (size(HOSB) == nSize)  return
         
       deallocate(HOSB)
   end if

   allocate(HOSB(nSize), stat=info)

   if (info/=0) allocHOSB=.FALSE.
end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function allocHOSBCX(nSize)
   integer, intent(IN) :: nSize

   integer :: info

   allocHOSBCX = .TRUE. 
 
   if (allocated(HOSBCX)) then
       if (size(HOSBCX) == nSize) return
       deallocate(HOSBCX)
   end if
   
   allocate(HOSBCX(nSize), stat=info)

   if (info/=0) allocHOSBCX=.FALSE.
end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function allocOSB()
    integer :: info

    call deallocOSB()

    allocOSB=.false.
    
    allocate(VOSB(myVOSB%mLen), RES(myLen), EIG0(myLen), stat=info)
    if (info/=0) return
    
    if ((sJOB==JOB_RES1) .OR.(sJOB==JOB_RES2) .OR. (sJOB==JOB_CRP) &
	.OR.(sJOB==JOB_CRP1) .OR. (sJOB==JOB_CRP2)  )then
        allocate(APP(myLen),APR(myLen),AP(myLen), stat=info)
    end if
    if (info/=0) return

    if (sCX) then
       allocate(H0CX(myH0%mLen), stat=info)
       if (info/=0) return
       
       if (sST) then
         allocate(HOSBCX(myHOSB%mLen), stat=info)
         if (info/=0) return
       end if

       if (sNDVR) then
          allocate(OUTHCX(outLen),stat=info)
          if (info/=0) return
       end if

       if (totalDep) then
          allocate(DEPCX(myDEP%mLen))
       end if
    else
       allocate(H0(myH0%mLen), stat=info)
       if (info/=0) return
       
       if (sST) then
         allocate(HOSB(myHOSB%mLen), stat=info)
         if (info/=0) return
       end if

       if (sNDVR) then
          allocate(OUTH(outLen),stat=info)
          if (info/=0) return
       end if

       if (totalDep) then
          allocate(DEP(myDEP%mLen))
       end if
    end if

    allocOSB = .true.
    
end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine deallocOSB()
   
    if (allocated(H0))     deallocate(H0)
    if (allocated(H0CX))   deallocate(H0CX)
    if (allocated(RES))    deallocate(RES)

    if (allocated(SQ_APR)) deallocate(SQ_APR)
    if (allocated(APP))    deallocate(APP)
    if (allocated(APR))    deallocate(APR)
    if (allocated(AP))     deallocate(AP)

    if (allocated(DEP))    deallocate(DEP)
    if (allocated(DEPCX))  deallocate(DEPCX)

    if (allocated(OUTH))   deallocate(OUTH)
    if (allocated(OUTHCX)) deallocate(OUTHCX)

    if (allocated(HOSB))   deallocate(HOSB)
    if (allocated(HOSBCX)) deallocate(HOSBCX)

    if (allocated(VOSB))   deallocate(VOSB)
    if (allocated(EIG0))   deallocate(EIG0)

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


