!ccccccccccccccccccccccccccccccccccccccccccccccccc 
!c      PX part of OSB_BASE:                     c
!c  This should be called after initOSBW() if    c
!c     osbw preconditioner is selected,          c
!ccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = V * ( H0)^-1 * V^T * X              c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine PX(NIN, X, Y)     
      integer, intent(IN)    :: NIN
      double precision, intent(IN)  :: X(NIN)
      double precision, intent(OUT) :: Y(NIN)

      call MYPX(sOSB, sPC, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = V * ( H0 - E )^-1 * V^T * X         c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine PX0(NIN, X, Y)
      integer, intent(IN)    :: NIN
      double precision, intent(IN) :: X(NIN)
      double precision, intent(OUT) :: Y(NIN)

      call MYPX(sOSB, TA0, NIN, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = V * ( H0 - E )^-1 * V^T * X         c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX(NIN, X, Y)
      integer, intent(IN)    :: NIN
      double precision, intent(IN) :: X(NIN)
      double precision, intent(OUT) :: Y(NIN)

      call MYPX(sOSB, TA1, NIN, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = V * ( H0 - E )^-1 * V^T * X         c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPX0(NIN, X, Y)
      integer, intent(IN)    :: NIN
      double precision, intent(IN) :: X(NIN)
      double precision, intent(OUT) :: Y(NIN)

      call MYPX(TOSB, TA1, NIN, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = V * ( H0 - E )^-1 * V^T * X               c
!c                1.0/(H0-E)   : H0 is not in energy window c
!c   (H0-E)^-1 =  -1.0/DE      : H0 is in [E-DELTA_E,E]     c
!c                 1.0/DE      : H0 is in [E, E+DELTA_E]    c
!c    X and Y should not be the same.                       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine EPXD1(NIN, X, Y)
      integer, intent(IN)    :: NIN
      double precision, intent(IN) :: X(NIN)
      double precision, intent(OUT) :: Y(NIN)

      call MYPX(TOSBD1,TA1, NIN, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = V * ( H0 - E )^-1 * V^T * X                  c 
!c                1.0/(H0-E)   : H0 is not in energy window    c
!c   (H0-E)^-1 =  -1.0/(sig(H0)*alpha + (1+beta)*H0)           c
!c                 H0 is in [E-DELTA_E, E+DELTA_E]             c
!c    X and Y should not be the same.                          c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPXD2(NIN, X, Y)
      integer, intent(IN)    :: NIN
      double precision, intent(IN) :: X(NIN)
      double precision, intent(OUT) :: Y(NIN)

      call MYPX(TOSBD2, TA1, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             Y = V * H1^-1 * V^T * X                     c
!c         X and Y should not be the same.                 c
!c     It should be called after createOSBWParam()         c    
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EPXW(NIN, X, Y)
      integer, intent(IN) :: NIN
      double precision, intent(IN)  :: X(NIN)
      double precision, intent(OUT) :: Y(NIN)
     
      call MYPX(TOSBW, TA1, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                  Do the real work                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine MYPX(nOSB, nPC, NIN, X, Y)
      integer, intent(IN) :: nOSB,nPC,NIN
      double precision, intent(IN)  :: X(NIN)
      double precision, intent(OUT) :: Y(NIN)

!***********************************************
      integer :: level, i, info, osbwcnt
      double precision :: DE, RDE, ADE, E0, tmp
      double precision :: alpha, beta

!**************************************

      Y(1:NIN) = X(1:NIN); 
      DE = DABS(sOSBW%mDE);   RDE  = 1.0D0 / DE

      do level = sF, 1, -1
          call VTX(myDim(level), sN(level), myBLK(level),   &
                   VOSB(myVOSB%mStart(level)),Y)
      end do

         ! print *, 'pX vector:',Y
         ! print *, 'Eig:',EIG0
      select case (nPC)
      case (TA0,TAAP0,TMAP0)
          E0 = 0.0D0
      case (TA1)
          E0 = sOSBW%mE0
      case default
          E0 = sOSBW%mE0
      end select

      PX_KERNEL: select case (nOSB)
      case (TOSB)     ! 1/Eig,  (Eig-E0)^-1
            Y(1:NIN) = Y(1:NIN)/( EIG0(1:NIN) - E0) 

      case (TOSBD1)         ! simple version of OSBD
           do I = 1, NIN
              ADE = EIG0(i) - E0
              if (DABS(ADE) <= DE) then
                 if (ADE<0.0D0) then 
                    Y(i) = -Y(i)*RDE
                 else
                    Y(i) = Y(i)*RDE
                 end if
              else          
                 Y(i) = Y(i)/ADE
              end if
           end do

      case (TOSBD2)      ! (alpha+(1+beta)*DE)^-1
           beta  = 1.0D0 + sOSBW%mBeta
           alpha = -sOSBW%mBeta*DE
           do I = 1, NIN
              ADE = EIG0(i) - E0
              if ( DABS(ADE) <= DE ) then
                  if (ADE<0.0D0) then
                      tmp = -alpha + beta*ADE
                  else
                      tmp = alpha + beta*ADE
                  end if
              else          
                  tmp = ADE
              end if
              Y(i) = Y(i)/tmp
           end do

      case (TOSBW)            
          do i = 1, hwLen 
             WY(i) = Y(WIND(i))
          end do
            ! print *, 'Initial vector:', WY 
          HWTMP(1:hwLen,1:hwLen) = HWR(1:hwLen,1:hwLen)
         
            ! call LAPACK subroutine to get A^-1 * X
          call DSYTRS('U',hwLen,1,HWTMP,hwLen,WIPIV,WY,hwLen,INFO) 
         
            !print *, 'final vector: A^-1*Y', WY 
          Y(1:NIN) = Y(1:NIN)/( EIG0(1:NIN) - E0)
            !print *, 'wind:', wind
          do i = 1, hwLen
             Y(WIND(i)) = WY(i)
          end do

      case default     ! 1/Eig,  (Eig-E0)^-1
            Y(1:NIN) = Y(1:NIN)/( EIG0(1:NIN) - E0)

      end select PX_KERNEL
 
       ! print *, 'After kernel: pX vector:',Y
      do level = 1, sF
         call VX(myDim(level),sN(level),myBLK(level),   &
                 VOSB(myVOSB%mStart(level)),Y)
      end do
      !print *, 'final pX vector:',Y
  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
