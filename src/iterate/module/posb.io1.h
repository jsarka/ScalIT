!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadOutH()
   integer :: layer

   layer=sF
   if (sCX) then
      loadOutH = loadGrid_CX(layer,fOutH,OutHCX)
   else
      loadOutH = loadGrid(layer,fOutH,OutH)
   end if
end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadRES()
   integer :: layer

   layer=1
   loadRES = loadSeq(layer,fRES,ResSeq)

   if (loadRES) then
      layer=sF
      loadRES = loadDiag(layer,fRES,RES)
   end if

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadEig()
   integer :: layer

   layer=1
   loadEig = loadSeq(layer,fEig,Eig0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveEig()
   integer :: layer

   layer=1
   saveEig = saveSeq(layer,fEig,Eig0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadAP()
   integer :: layer

   layer=1
   loadAP = loadSeq(layer,fAPP,AP)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadAPP()
   integer :: layer

   layer=1
   loadAPP = loadSeq(layer,fAPP,APP)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadAPR()
   integer :: layer

   layer=1
   loadAPR = loadSeq(layer,fAPR,APR)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveVOSBEig()

     saveVOSBEig = saveVOSB()

     if (saveVOSBEig) saveVOSBEig = saveEig()

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadVOSBEig()

     loadVOSBEig = loadVOSB()

     if (loadVOSBEig) loadVOSBEig = loadEig()

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      Not implemented
logical function loadDep()
   integer :: layer

   layer=1
   loadDep = .true.

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c The following subroutines are not implemented without parallel IO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadHOSB()
    integer :: i
  
    loadHOSB = .true.

    if (.NOT. sST) return

    if (sCX) then

    else

    end if

end function
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveHOSB()
 
    saveHOSB=.true.

    if (.NOT. sST) return

end function
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadVOSB()


    loadVOSB = .true.
end function
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveVOSB()

    saveVOSB = .true.
end function
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
