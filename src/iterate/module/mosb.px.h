!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!c     Preconditioner Matrix-vector product: Real version       c
!c     This should be called after initOSBW() if osbw           c
!c   preconditioner is selected.  X, Y should not be the same.  c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = V * ( H0)^-1 * V^T * X              c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine PX(N, X, Y)     
      integer, intent(IN)    :: N
      double precision, intent(IN)  :: X(N)
      double precision, intent(OUT) :: Y(N)

      call MYPX(sOSB, sPC, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = V * ( H0 - E )^-1 * V^T * X         c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine PX0(N, X, Y)
      integer, intent(IN)    :: N
      double precision, intent(IN) :: X(N)
      double precision, intent(OUT) :: Y(N)

      call MYPX(sOSB, TA0, N, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = V * ( H0 - E )^-1 * V^T * X         c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX(N, X, Y)
      integer, intent(IN)    :: N
      double precision, intent(IN) :: X(N)
      double precision, intent(OUT) :: Y(N)

      call MYPX(sOSB, TA1, N, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = V * ( H0 - E )^-1 * V^T * X         c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX0(N, X, Y)
      integer, intent(IN)    :: N
      double precision, intent(IN) :: X(N)
      double precision, intent(OUT) :: Y(N)

      call MYPX(TOSB, TA1, N, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = V * ( H0 - E )^-1 * V^T * X               c
!c                1.0/(H0-E)   : H0 is not in energy window c
!c   (H0-E)^-1 =  -1.0/DE      : H0 is in [E-DELTA_E,E]     c
!c                 1.0/DE      : H0 is in [E, E+DELTA_E]    c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine EPXD1(N, X, Y)
      integer, intent(IN)    :: N
      double precision, intent(IN) :: X(N)
      double precision, intent(OUT) :: Y(N)

      call MYPX(TOSBD1,TA1, N, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = V * ( H0 - E )^-1 * V^T * X                  c 
!c                1.0/(H0-E)   : H0 is not in energy window    c
!c   (H0-E)^-1 =  -1.0/(sig(H0)*alpha + (1+beta)*H0)           c
!c                 H0 is in [E-DELTA_E, E+DELTA_E]             c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPXD2(N, X, Y)
      integer, intent(IN)    :: N
      double precision, intent(IN) :: X(N)
      double precision, intent(OUT) :: Y(N)

      call MYPX(TOSBD2, TA1, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             Y = V * H1^-1 * V^T * X                     c
!c     It should be called after createOSBWParam()         c    
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPXW(N, X, Y)
      integer, intent(IN) :: N
      double precision, intent(IN)  :: X(N)
      double precision, intent(OUT) :: Y(N)
     
      call MYPX(TOSBW, TA1, N, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
