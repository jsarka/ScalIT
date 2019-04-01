!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                     Initialize HOSB data                            c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine myHOSBInit(level)
   integer, intent(IN) :: level

   call HOSBINIT_XYZ(myBlk(level), sN(level), myDim(level),       &
         H0(myH0%mStart(level)), EIG0, HOSB(myHOSB%mStart(level)))

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine myHOSBInit_CX(level)
   integer, intent(IN) :: level

   call HOSBINIT_XYZ_CX(myBlk(level), sN(level), myDim(level),     &
        H0CX(myH0%mStart(level)), EIG0, HOSBCX(myHOSB%mStart(level)))           

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine myHOSBDepInit(level)
   integer, intent(IN) :: level

   call HOSBINIT_DEP(myBlk(level),sN(level),myDim(level),        & 
            H0(myH0%mStart(level)),DEP(myDEP%mStart(level)),     &
              EIG0, HOSB(myHOSB%mStart(level)))

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine myHOSBDepInit_CX(level)
   integer, intent(IN) :: level
        
   call HOSBINIT_DEP_CX(myBlk(level),sN(level),myDim(level),       & 
             H0CX(myH0%mStart(level)), DEPCX(myDEP%mStart(level)), &
               EIG0, HOSBCX(myHOSB%mStart(level)))

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

