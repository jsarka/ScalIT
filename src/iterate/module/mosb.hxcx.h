!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Original Matrix-Vector product:HX              c
!c  H: Complex, X,Y: complex, X, Y should be different  c
!c           Y = H*X, (H-E)*X, (H(+/-)iAp-E)*X          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  HX_CX(N, X, Y)
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

      call MYHX_CX(sHC, N, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  HijX_CX(N, X, Y)
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

!      call MYHX_CX(sHij, N, X, Y)
      call MYHX_CX(TA0, N, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y =  H  * X                             c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  HX0_CX(N, X, Y)
      integer, intent(IN)    :: N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

      call MYHX_CX(TA0, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H - E ) * X                       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_CX(N, X, Y)    
      integer, intent(IN)    :: N     
      double complex,intent(IN)  :: X(N)  
      double complex,intent(OUT) :: Y(N)

      call MYHX_CX(TA1, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y =  ( H + iAP) * X                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  AHX_CX(N, X, Y)
      integer, intent(IN)    :: N 
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)
      
      call MYHX_CX(TAAP0, N, X, Y)
  
  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y =  ( H - iAP )  * X                   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  MHX_CX(N, X, Y)
      integer, intent(IN)    :: N 
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)
      
      call MYHX_CX(TMAP0, N, X, Y)
  
  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H + iAP - E ) * X                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_PCX(N, X, Y)       
      integer, intent(IN)    :: N     
      double complex,intent(IN)  :: X(N)   
      double complex,intent(OUT) :: Y(N)   

      call MYHX_CX(TAAP, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H + iAPP - E ) * X                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_PPCX(N, X, Y)       
      integer, intent(IN)    :: N     
      double complex,intent(IN)  :: X(N) 
      double complex,intent(OUT) :: Y(N) 

      call MYHX_CX(TAAPP, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H + iAPR - E ) * X                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_PRCX(N, X, Y)       
      integer, intent(IN)    :: N     
      double complex,intent(IN)  :: X(N) 
      double complex,intent(OUT) :: Y(N) 

      call MYHX_CX(TAAPR, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H - iAP - E ) * X                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_MCX(N, X, Y)    
      integer, intent(IN)    :: N     
      double complex,intent(IN)  :: X(N)
      double complex,intent(OUT) :: Y(N)

      call MYHX_CX(TMAP, N, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H - iAPP - E ) * X                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_MPCX(N, X, Y)    
      integer, intent(IN)    :: N     
      double complex,intent(IN)  :: X(N)
      double complex,intent(OUT) :: Y(N)

      call MYHX_CX(TMAPP, N, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H - iAPR - E ) * X                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_MRCX(N, X, Y)    
      integer, intent(IN)    :: N     
      double complex,intent(IN)  :: X(N)
      double complex,intent(OUT) :: Y(N)

      call MYHX_CX(TMAPR, N, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
