!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!c   Preconditioner Matrix-vector product: Complex version    c
!c     This should be called after initOSBW() if osbw         c
!c    preconditioner is selected.X,Y should not be the same.  c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = V * ( H0)^-1 * V^T * X              c
!c          X and Y should not be the same.           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine PX_DX(N, X, Y)
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

      call MYPX_DX(sOSB, sPC, N, X, Y)
 
  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine PX0_DX(N, X, Y)           
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

      call MYPX_DX(sOSB,TA0, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c          Y = V * ( H0 - E )^-1 * V^T * X           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_DX(N, X, Y)              ! OSB
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

      call MYPX_DX(sOSB,TA1, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Y = V * ( H0 - E +iAP )^-1 * V^T * X         c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_PCX(N, X, Y)              ! OSB
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

      call MYPX_DX(sOSB,TAAP, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Y = V * ( H0 - E -iAP )^-1 * V^T * X           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_MCX(N, X, Y)              ! OSB
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

      call MYPX_DX(sOSB,TMAP, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c          Y = V * ( H0 - E )^-1 * V^T * X           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_DX0(N, X, Y)              ! OSB
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

      call MYPX_DX(TOSB,TA1, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_DXD1(N, X, Y)            ! OSBD1
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

      call MYPX_DX(TOSBD1,TA1, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_DXD2(N, X, Y)            ! OSBD2
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

      call MYPX_DX(TOSBD2,TA1, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_DXW(N, X, Y)             ! OSBW
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

      call MYPX_DX(TOSBW, TA1, N, X, Y)

   end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Y = V * ( H0 - E + iAP )^-1 * V^T * X        c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_PCX0(N, X, Y)             ! OSB
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

      call MYPX_DX(TOSB,TAAP, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_PCXD1(N, X, Y)           ! OSBD1
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

      call MYPX_DX(TOSBD1,TAAP, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_PCXD2(N, X, Y)           ! OSBD2
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

      call MYPX_DX(TOSBD2,TAAP, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_PCXW(N, X, Y)           ! OSBW
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

      call MYPX_DX(TOSBW, TAAP, N, X, Y)

   end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Y = V * ( H0 - E - iAP)^-1 * V^T * X         c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_MCX0(N, X, Y)             ! OSB
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

      call MYPX_DX(TOSB,TMAP, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_MCXD1(N, X, Y)           ! OSBD1
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

      call MYPX_DX(TOSBD1,TMAP, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_MCXD2(N, X, Y)           ! OSBD2
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

      call MYPX_DX(TOSBD2,TMAP, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX_MCXW(N, X, Y)           ! OSBW
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

      call MYPX_DX(TOSBW, TMAP, N, X, Y)

   end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
