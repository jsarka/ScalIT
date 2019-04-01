!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  HX part of OSB     :                             c
!c  Subroutines:                                     c
!c  HX    (NIN, X, Y):  Y = H * X                    c
!c  EHX   (NIN, X, Y):  Y = ( H - E) * X             c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  HX_DX(NIN, X, Y)
      integer, intent(IN)    :: NIN
      double complex, intent(IN)  :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      call MYHX_DX(sHC, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccc
!c             Y = H * X                    c
!c        X, Y should be different          c
!cccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  HX0_DX(NIN, X, Y)
      integer, intent(IN)    :: NIN
      double complex, intent(IN)  :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      call MYHX_DX(TA0, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H - E ) * X                       c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_DX(NIN, X, Y)    
      integer, intent(IN)    :: NIN     
      double complex,intent(IN)  :: X(NIN)  
      double complex,intent(OUT) :: Y(NIN)

      call MYHX_DX(TA1, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H + iAP ) * X                     c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  AHX_DX(NIN, X, Y)    
      integer, intent(IN)    :: NIN     
      double complex,intent(IN)  :: X(NIN)  
      double complex,intent(OUT) :: Y(NIN)

      call MYHX_DX(TAAP0, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H - iAP ) * X                     c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  MHX_DX(NIN, X, Y)    
      integer, intent(IN)    :: NIN     
      double complex,intent(IN)  :: X(NIN)  
      double complex,intent(OUT) :: Y(NIN)

      call MYHX_DX(TMAP0, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H + iAP - E ) * X                 c
!c    This function is the same as EHX_CX             c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_PDX(NIN, X, Y)       
      integer, intent(IN)    :: NIN     
      double complex,intent(IN)  :: X(NIN)   
      double complex,intent(OUT) :: Y(NIN)

      call MYHX_DX(TAAP, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H + iAPP - E ) * X                c
!c    This function is the same as EHX_DX             c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_PPDX(NIN, X, Y)       
      integer, intent(IN)    :: NIN     
      double complex,intent(IN)  :: X(NIN)   
      double complex,intent(OUT) :: Y(NIN)

      call MYHX_DX(TAAPP, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H + iAPR - E ) * X                 c
!c    This function is the same as EHX_DX             c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_PRDX(NIN, X, Y)       
      integer, intent(IN)    :: NIN     
      double complex,intent(IN)  :: X(NIN)   
      double complex,intent(OUT) :: Y(NIN)

      call MYHX_DX(TAAPR, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H - iAP - E ) * X                 c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_MDX(NIN, X, Y)    
      integer, intent(IN)    :: NIN     
      double complex,intent(IN)  :: X(NIN)
      double complex,intent(OUT) :: Y(NIN)

      call MYHX_DX(TMAP, NIN, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H - iAPP - E ) * X                c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_MPDX(NIN, X, Y)    
      integer, intent(IN)    :: NIN     
      double complex,intent(IN)  :: X(NIN)
      double complex,intent(OUT) :: Y(NIN)

      call MYHX_DX(TMAPP, NIN, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H - iAPR - E ) * X                c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_MRDX(NIN, X, Y)    
      integer, intent(IN)    :: NIN     
      double complex,intent(IN)  :: X(NIN)
      double complex,intent(OUT) :: Y(NIN)

      call MYHX_DX(TMAPR, NIN, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Do the real work                        c
!c               Y = (H-E)*X or H*X                   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  MYHX_DX(nType, NIN, X, Y)    
      integer, intent(IN)    :: nType, NIN
      double complex,intent(IN)  :: X(NIN) 
      double complex,intent(OUT) :: Y(NIN)

      double complex :: tmp(NIN)
      double precision :: E0
      integer :: level

      E0 = sOSBW%mE0

      select case (nType)
      case (TA0)    !  H*X
           Y(1:NIN) = RES(1:NIN) * X (1:NIN)

      case (TA1)    !  (H-E)*X
           Y(1:NIN) = (RES(1:NIN)-E0) * X (1:NIN)

      case (TAAP0)  !  (H+iAP)*X
           Y(1:NIN) = DCMPLX(RES(1:NIN),AP(1:NIN)) * X (1:NIN)

      case (TMAP0)  !  (H-iAP)*X 
           Y(1:NIN) = DCMPLX(RES(1:NIN),-AP(1:NIN)) * X (1:NIN)

      case (TAAP)   !  (H+iAP-E)*X
           Y(1:NIN) = DCMPLX(RES(1:NIN)-E0,AP(1:NIN)) * X (1:NIN)

      case (TAAPP)  !  (H+iAPP-E)*X
           Y(1:NIN) = DCMPLX(RES(1:NIN)-E0,APP(1:NIN)) * X (1:NIN)

      case (TAAPR)  !  (H+iAPR-E)*X
           Y(1:NIN) = DCMPLX(RES(1:NIN)-E0,APR(1:NIN)) * X (1:NIN)

      case (TMAP)   !  (H-iAP-E)*X         
           Y(1:NIN) = DCMPLX(RES(1:NIN)-E0,-AP(1:NIN)) * X (1:NIN)

      case (TMAPP)  !  (H-iAPP-E)*X
           Y(1:NIN) = DCMPLX(RES(1:NIN)-E0,-APP(1:NIN)) * X (1:NIN)

      case (TMAPR)  !  (H-iAPR-E)*X
           Y(1:NIN) = DCMPLX(RES(1:NIN)-E0,-APR(1:NIN)) * X (1:NIN)

      case default   !default, (H-E)*X
           Y(1:NIN) = (RES(1:NIN)-E0) * X (1:NIN)

      end select

      do level = 1, sF
         if (sNDVR .and. (level==sF)) then       
  	    call H3X_DX(myDim(level),sN(level),OUTH, X, tmp)
         else
      	    if (sDEP(level)) then
                call H1X_DEP_DX (myDim(level), sN(level), myBLK(level),     &
                     H0(myH0%mStart(level)),DEP(myDEP%mStart(level)), X, tmp)
            else
                call H1X_XYZ_DX (myDim(level), sN(level), myBLK(level),     &
                        H0(myH0%mStart(level)), X, tmp)
            end if
         end if

         Y(1:NIN) = Y(1:NIN) + tmp(1:NIN)     
       end do
  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
