!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Original Matrix-Vector product:HX              c
!c     H: Real, X,Y: Real, X, Y should be different     c
!c                   Y = H*X, (H-E)*X                   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine  HX(N, X, Y)
      integer, intent(IN)    :: N
      double precision,intent(IN)  :: X(N)
      double precision,intent(OUT) :: Y(N)

      call MYHX(sHC, N, X, Y)

 end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine  HijX(N, X, Y)
      integer, intent(IN)    :: N
      double precision,intent(IN)  :: X(N)
      double precision,intent(OUT) :: Y(N)

!      call MYHX(sHij, N, X, Y)
      call MYHX(TA0, N, X, Y)

 end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H - E ) * X                       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine  HX0(N, X, Y)
      integer, intent(IN)    :: N
      double precision,intent(IN)  :: X(N)
      double precision,intent(OUT) :: Y(N)

      call MYHX(TA0, N, X, Y)

 end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H - E ) * X                       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX(N, X, Y)    
      integer, intent(IN)    :: N
      double precision,intent(IN)  :: X(N)
      double precision,intent(OUT) :: Y(N)

      call MYHX(TA1, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



