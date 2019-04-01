!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Original Matrix-Vector product:HX              c
!c    H: Real, X,Y: complex, X, Y should be different   c
!c           Y = H*X, (H-E)*X, (H(+/-)iAp-E)*X          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  HX_DX(N, X, Y)
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

      call MYHX_DX(sHC, N, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  HijX_DX(N, X, Y)
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

!      call MYHX_DX(sHij, N, X, Y)
      call MYHX_DX(TA0, N, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c               Y =  H  * X                          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  HX0_DX(N, X, Y)
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

      call MYHX_DX(TA0, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H - E ) * X                       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_DX(N, X, Y)    
      integer, intent(IN)    :: N     
      double complex,intent(IN)  :: X(N)  
      double complex,intent(OUT) :: Y(N)

      call MYHX_DX(TA1, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H + iAP ) * X                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  AHX_DX(N, X, Y)    
      integer, intent(IN)    :: N     
      double complex,intent(IN)  :: X(N)
      double complex,intent(OUT) :: Y(N)
      
      call MYHX_DX(TAAP0, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H - iAP ) * X                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  MHX_DX(N, X, Y)    
      integer, intent(IN)    :: N     
      double complex,intent(IN)  :: X(N)
      double complex,intent(OUT) :: Y(N)
      
      call MYHX_DX(TMAP0, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H + iAP - E ) * X                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_PDX(N, X, Y)       
      integer, intent(IN)    :: N     
      double complex,intent(IN)  :: X(N)   
      double complex,intent(OUT) :: Y(N)   

      call MYHX_DX(TAAP, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H + iAPP - E ) * X                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_PPDX(N, X, Y)       
      integer, intent(IN)    :: N     
      double complex,intent(IN)  :: X(N)   
      double complex,intent(OUT) :: Y(N)   

      call MYHX_DX(TAAPP, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H + iAPR - E ) * X                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_PRDX(N, X, Y)       
      integer, intent(IN)    :: N     
      double complex,intent(IN)  :: X(N)   
      double complex,intent(OUT) :: Y(N)   

      call MYHX_DX(TAAPR, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H - iAP - E ) * X                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_MDX(N, X, Y)    
      integer, intent(IN)    :: N     
      double complex,intent(IN)  :: X(N)
      double complex,intent(OUT) :: Y(N)

      call MYHX_DX(TMAP, N, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H - iAPP - E ) * X                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_MPDX(N, X, Y)    
      integer, intent(IN)    :: N     
      double complex,intent(IN)  :: X(N)
      double complex,intent(OUT) :: Y(N)

      call MYHX_DX(TMAPP, N, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H - iAPR - E ) * X                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_MRDX(N, X, Y)    
      integer, intent(IN)    :: N     
      double complex,intent(IN)  :: X(N)
      double complex,intent(OUT) :: Y(N)

      call MYHX_DX(TMAPR, N, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



