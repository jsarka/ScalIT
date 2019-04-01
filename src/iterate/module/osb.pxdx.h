!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!c  PX part of OSB_BASE:Complex version                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = V * ( H0)^-1 * V^T * X              c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine PX_DX(NIN, X, Y)
      integer, intent(IN)    :: NIN
      double complex, intent(IN) :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      call MYPX_DX(sOSB, sPC, NIN, X, Y)
 
  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = V * ( H0 - E )^-1 * V^T * X         c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine PX0_DX(NIN, X, Y)           ! OSB
      integer, intent(IN)    :: NIN
      double complex, intent(IN) :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      call MYPX_DX(sOSB,TA0, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_DX(NIN, X, Y)           ! OSB
      integer, intent(IN)    :: NIN
      double complex, intent(IN) :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      call MYPX_DX(sOSB,TA1, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         Y = V * ( H0 - E + iAP )^-1 * V^T * X      c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_PCX(NIN, X, Y)
      integer, intent(IN)    :: NIN
      double complex, intent(IN) :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      call MYPX_DX(sOSB, TAAP, NIN, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         Y = V * ( H0 - E - iAP)^-1 * V^T * X       c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_MCX(NIN, X, Y)
      integer, intent(IN)    :: NIN
      double complex, intent(IN) :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      call MYPX_DX(sOSB, TMAP, NIN, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_DX0(NIN, X, Y)           ! OSB
      integer, intent(IN)    :: NIN
      double complex, intent(IN) :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      call MYPX_DX(TOSB,TA1, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_DXD1(NIN, X, Y)           ! OSBD1
      integer, intent(IN)    :: NIN
      double complex, intent(IN) :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      call MYPX_DX(TOSBD1,TA1, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_DXD2(NIN, X, Y)           ! OSBD2
      integer, intent(IN)    :: NIN
      double complex, intent(IN) :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      call MYPX_DX(TOSBD1,TA1, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c              Full  version  of  OSBW                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_DXW(NIN, X, Y)
      integer, intent(IN)    :: NIN
      double complex, intent(IN) :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      call MYPX_DX(TOSBW, TA1, NIN, X, Y)

   end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_PCX0(NIN, X, Y)           ! OSB
      integer, intent(IN)    :: NIN
      double complex, intent(IN) :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      call MYPX_DX(TOSB,TAAP, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_PCXD1(NIN, X, Y)           ! OSBD1
      integer, intent(IN)    :: NIN
      double complex, intent(IN) :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      call MYPX_DX(TOSBD1,TAAP, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_PCXD2(NIN, X, Y)           ! OSBD2
      integer, intent(IN)    :: NIN
      double complex, intent(IN) :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      call MYPX_DX(TOSBD1,TAAP, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c              Full  version  of  OSBW                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_PCXW(NIN, X, Y)
      integer, intent(IN)    :: NIN
      double complex, intent(IN) :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      call MYPX_DX(TOSBW, TAAP, NIN, X, Y)

   end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_MCX0(NIN, X, Y)           ! OSB
      integer, intent(IN)    :: NIN
      double complex, intent(IN) :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      call MYPX_DX(TOSB,TMAP, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_MCXD1(NIN, X, Y)           ! OSBD1
      integer, intent(IN)    :: NIN
      double complex, intent(IN) :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      call MYPX_DX(TOSBD1,TMAP, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_MCXD2(NIN, X, Y)           ! OSBD2
      integer, intent(IN)    :: NIN
      double complex, intent(IN) :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      call MYPX_DX(TOSBD1,TMAP, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c              Full  version  of  OSBW                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_MCXW(NIN, X, Y)
      integer, intent(IN)    :: NIN
      double complex, intent(IN) :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      call MYPX_DX(TOSBW, TMAP, NIN, X, Y)

   end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
