!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                     Initialize HOSB data                            c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine myHOSBInit(level)
   integer, intent(IN) :: level

   call HOSBINIT_XYZ(blk(level), nout(level), nin(level),          &
         H0(myH0%mStart(level)), EIG0, HOSB(myHOSB%pStart(level)))

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine myHOSBInit_CX(level)
   integer, intent(IN) :: level

   call HOSBINIT_XYZ_CX(blk(level), nout(level), nin(level),       &
        H0CX(myH0%mStart(level)), EIG0, HOSBCX(myHOSB%pStart(level)))           

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine myHOSBDepInit(level)
   integer, intent(IN) :: level

   call HOSBINIT_DEP(blk(level),nout(level),nin(level),          & 
            H0(myH0%mStart(level)),DEP(myDEP%pStart(level)),     &
            EIG0, HOSB(myHOSB%pStart(level)))

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine myHOSBDepInit_CX(level)
   integer, intent(IN) :: level
        
   call HOSBINIT_DEP_CX(blk(level),nout(level),nin(level),         & 
           H0CX(myH0%mStart(level)), DEPCX(myDEP%pStart(level)),   &
           EIG0, HOSBCX(myHOSB%pStart(level)))

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

