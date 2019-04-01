!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function initMOSBW()
    initMOSBW = preInitMOSBW()


    if (initMOSBW) initMOSBW = postInitMOSBW()

    !if (id==2) then
    !  print *, 'OSBW:len=',hwLen
    !  print *,'HW:', HW 
    !  print *, 'HWR:', HWR
    !endif

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine finalMOSBW()
   
    call DeallocHW()

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function preInitMOSBW()

   integer :: i
   call adjustWinSize(sOSBW%mCnt, sOSBW%mE0, sOSBW%mDE)

   call getWGlobalSize(sOSBW%mE0, sOSBW%mDE, phwLen, hwLen)

   if (id==rootID) print *, 'hwLen', hwLen

   preInitMOSBW = .false.

   if (hwlen==0) then
      sPC = sPC - 2;  return      ! use OSBD
   else
      preInitMOSBW = allocHW()
   end if

   if (preInitMOSBW) then
       if (phwLen>0) then 
           if (id==rootID) then
               print *, "    Entering getWGlobalIndex with phwLen = ", phwLen  
           end if
           call getWGlobalIndex(sOSBW%mE0,sOSBW%mDE,phwLen,hwLen,gCnt,  &
                          sCnt,pInd,gInd,nodInd)
       else
           call getWGlobalIndex(sOSBW%mE0,sOSBW%mDE, 1,   hwLen,gCnt,   &
                          sCnt,pInd,gInd,nodInd)
       end if

       call calGlobalIndex(hwLen, gInd, gKIndex, gPIndex)

       call getWinVector(Eig0, gEig0)

       if (sAP) then
          call getWinVector(APP,gAPP)
          call getWinVector(APR,gAPR)
          gAP(1:hwLen) = gAPP(1:hwLen) + gAPR(1:hwLen)
       end if

       do i = 1, hwLen
          call mgetViColIndex(nNodes,sF,sN,myconf%gDim,myconf%gBlk, &
                     gKIndex(i), grpInd(1,i), rootInd(1,i),         &
                     nodeNum(1,i), blkInd(1,i), sNInd(1,i) )
       end do

   end if
end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function postInitMOSBW()
 
    if (sHW>=0) then
       if (sCX) then
          call getMOSBWHCX(hwLen,gKIndex,gPIndex,HWCX)
       else
          call getMOSBWH(hwLen,gKIndex,gPIndex,HW)
       end if
       if (sHW>0) call saveHW()  
    else
       call loadHW()
    end if

    if (sCX) then
        postInitMOSBW = createWHCX()        
    else
        if (sAP) then
           postInitMOSBW = createWHDX()
        else
           postInitMOSBW = createWH()
        end if
    end if

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function allocHW()

   integer :: info, tmpLen
  
   allocHW = .false.

   if (hwLen<=0) return

   call deallocHW()

   if (sCX) then
      allocate(HWCX(hwLen,hwLen), HWRCX(hwLen,hwLen),        &
               HWTMPCX(hwLen,hwLen),WXCX(hwLen),WYCX(hwLen), &
               stat=info)
   else
      allocate(HW(hwLen,hwLen),HWR(hwLen,hwLen),             &
              HWTMP(hwLen,hwLen),WX(hwLen),WY(hwLen),stat=info)
   end if
   if (info /=0) return
 
   allocate(WIPIV(hwLen), WIND(hwLen), stat=info)
   if (info /= 0) return

   if (phwLen<=0) then
      tmpLen = 1
   else
      tmpLen = phwLen
   end if

   allocate(gCnt(nNodes),sCnt(nNodes),pInd(tmpLen), gInd(hwLen),   &
            gKIndex(hwLen),gPIndex(hwlen),nodInd(hwLen),stat=info)
   if (info/=0) return

   allocate(gAP(hwLen),gAPP(hwLen),gAPR(hwLen),gEig0(hwLen),stat=info)
   if (info/=0) return

   allocate(grpInd(sF,hwLen),rootInd(sF,hwLen),nodeNum(sF,hwLen),  &
            blkInd(sF,hwLen), sNInd(sF,hwLen), stat=info)

   allocHW = (info==0)

end function
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine deallocHW()

    if (allocated(HW))     deallocate(HW)
    if (allocated(HWR))    deallocate(HWR)
    if (allocated(HWTMP))  deallocate(HWTMP)

    if (allocated(WY))     deallocate(WY)
    if (allocated(WYCX))   deallocate(WYCX)

    if (allocated(HWCX))   deallocate(HWCX)
    if (allocated(HWRCX))  deallocate(HWRCX)
    if (allocated(HWTMPCX))   deallocate(HWTMPCX)

    if (allocated(WX))     deallocate(WX)
    if (allocated(WXCX))   deallocate(WXCX)

    if (allocated(WIPIV))  deallocate(WIPIV)
    if (allocated(WIND))   deallocate(WIND)

    if (allocated(gCnt))  deallocate(gCnt)
    if (allocated(sCnt))  deallocate(sCnt)

    if (allocated(pInd))  deallocate(pInd)
    if (allocated(gInd))  deallocate(gInd)

    if (allocated(gKIndex)) deallocate(gKIndex)
    if (allocated(gPIndex)) deallocate(gPIndex)
    if (allocated(nodInd))  deallocate(nodInd)

    if (allocated(gAP))     deallocate(gAP)
    if (allocated(gAPP))    deallocate(gAPP)
    if (allocated(gAPR))    deallocate(gAPR)
    if (allocated(gEig0))   deallocate(gEig0)

    if (allocated(grpInd))     deallocate(grpInd)
    if (allocated(rootInd))    deallocate(rootInd)
    if (allocated(nodeNum))    deallocate(nodeNum)
    if (allocated(blkInd))     deallocate(blkInd)
    if (allocated(sNInd))      deallocate(sNInd)

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

