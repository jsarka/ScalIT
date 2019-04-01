!ccccccccccccccccccccccccccccccccccccccccccccccc 
!c           VX part of OSB:                   c
!c   Calculate V*X or V^T*X at one level       c
!c  Subroutines:                               c
!c  VOSBX    (level, X)    !  X = V * X        c
!c  VOSBTX   (level, X)    !  X = V^T * X      c
!c  VOSBX_CX (level, X)    !  X = V * X        c
!c  VOSBTX_CX(level, X)    !  X = V^T * X      c
!ccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         VOSBX: Do the real V*X job at specific level        c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine VOSBX(level, X)
      integer, intent(IN) :: level
      double precision,intent(INOUT) :: X(myLen)

      call VX(myDim(level),sN(level),myBLK(level),     &
              VOSB(myVOSB%mStart(level)),X)
  
  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     VOSBTX: Do the real V^T*X job at specific level       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine VOSBTX(level, X)
      integer, intent(IN) :: level
      double precision,intent(INOUT) :: X(myLen)
 
      call VTX(myDim(level),sN(level),myBLK(level),   &
               VOSB(myVOSB%mStart(level)),X)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         V_MULT_X: Do the real V*X job at specific level        c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine VOSBX_CX(level, X)
      integer, intent(IN) :: level
      double complex, intent(INOUT) :: X(myLen)

      call VX_DX(myDim(level),sN(level),myBLK(level),     &
                 VOSB(myVOSB%mStart(level)),X)
  
  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     VT_MULT_X: Do the real V^T*X job at specific level       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine VOSBTX_CX(level, X)
      integer, intent(IN) :: level
      double complex, intent(INOUT) :: X(myLen)
 
      call VTX_DX(myDim(level),sN(level),myBLK(level),    &
                  VOSB(myVOSB%mStart(level)),X)

 end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
