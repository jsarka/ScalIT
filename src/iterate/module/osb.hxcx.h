!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  HX part of OSB     :                             c
!c  HX    (NIN, X, Y):  Y = H * X                    c
!c  EHX   (NIN, X, Y):  Y = ( H - E) * X             c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccc
!c             Y = H * X                    c
!c        X, Y should be different          c
!cccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  HX_CX(NIN, X, Y)
      integer, intent(IN)    :: NIN
      double complex, intent(IN)  :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      call MYHX_CX(sHC, NIN, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  HijX_CX(NIN, X, Y)
      integer, intent(IN)    :: NIN
      double complex, intent(IN)  :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      !call MYHX_CX(sHij, NIN, X, Y)
      call MYHX_CX(TA0, NIN, X, Y)
  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  HX0_CX(NIN, X, Y)
      integer, intent(IN)    :: NIN
      double complex, intent(IN)  :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

      call MYHX_CX(TA0, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H - E ) * X                       c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_CX(NIN, X, Y)    
      integer, intent(IN)    :: NIN     
      double complex,intent(IN)  :: X(NIN)  
      double complex,intent(OUT) :: Y(NIN)

      call MYHX_CX(TA1, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c              Y = ( H + iAP ) * X                   c
!c    This function is the same as EHX_CX             c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  AHX_CX(NIN, X, Y)       
      integer, intent(IN)    :: NIN     
      double complex,intent(IN)  :: X(NIN)   
      double complex,intent(OUT) :: Y(NIN)

      call MYHX_CX(TAAP0, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c              Y = ( H - iAP ) * X                   c
!c    This function is the same as EHX_CX             c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  MHX_CX(NIN, X, Y)       
      integer, intent(IN)    :: NIN     
      double complex,intent(IN)  :: X(NIN)   
      double complex,intent(OUT) :: Y(NIN)

      call MYHX_CX(TMAP0, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H + iAP - E ) * X                 c
!c    This function is the same as EHX_CX             c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_PCX(NIN, X, Y)       
      integer, intent(IN)    :: NIN     
      double complex,intent(IN)  :: X(NIN)   
      double complex,intent(OUT) :: Y(NIN)

      call MYHX_CX(TAAP, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H + iAPP - E ) * X                c
!c    This function is the same as EHX_CX             c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_PPCX(NIN, X, Y)       
      integer, intent(IN)    :: NIN     
      double complex,intent(IN)  :: X(NIN)   
      double complex,intent(OUT) :: Y(NIN)

      call MYHX_CX(TAAPP, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H + iAPR - E ) * X                c
!c    This function is the same as EHX_CX             c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_PRCX(NIN, X, Y)       
      integer, intent(IN)    :: NIN     
      double complex,intent(IN)  :: X(NIN)   
      double complex,intent(OUT) :: Y(NIN)

      call MYHX_CX(TAAPR, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H - iAP - E ) * X                 c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_MCX(NIN, X, Y)    
      integer, intent(IN)    :: NIN     
      double complex,intent(IN)  :: X(NIN)
      double complex,intent(OUT) :: Y(NIN)

      call MYHX_CX(TMAP, NIN, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H - iAPP - E ) * X                c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_MPCX(NIN, X, Y)    
      integer, intent(IN)    :: NIN     
      double complex,intent(IN)  :: X(NIN)
      double complex,intent(OUT) :: Y(NIN)

      call MYHX_CX(TMAPP, NIN, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H - iAPR - E ) * X                c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX_MRCX(NIN, X, Y)    
      integer, intent(IN)    :: NIN     
      double complex,intent(IN)  :: X(NIN)
      double complex,intent(OUT) :: Y(NIN)

      call MYHX_CX(TMAPR, NIN, X, Y)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Do the real work                        c
!c               Y = (H-E)*X or H*X                   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  MYHX_CX(nType, NIN, X, Y)    
      integer, intent(IN)    :: nType, NIN
      double complex,intent(IN)  :: X(NIN) 
      double complex,intent(OUT) :: Y(NIN)

      double complex   :: tmp(NIN)
      double precision :: E0
      integer :: level

      E0 = sOSBW%mE0

      select case (nType)

      case (TA0)     !  H*X        
            Y(1:NIN) = RES(1:NIN) * X (1:NIN)

      case (TA1)     !  (H-E)*X        
            Y(1:NIN) = (RES(1:NIN)-E0) * X (1:NIN)

      case (TAAP0)   !  (H+iAP)*X        
            Y(1:NIN) = DCMPLX(RES(1:NIN),AP(1:NIN)) * X (1:NIN)

      case (TMAP0)   !  (H-iAP)*X        
            Y(1:NIN) = DCMPLX(RES(1:NIN),-AP(1:NIN)) * X (1:NIN)

      case (TAAP)    !  (H+iAP-E)*X        
            Y(1:NIN) = DCMPLX(RES(1:NIN)-E0,AP(1:NIN)) * X (1:NIN)

      case (TAAPP)   !  (H+iAPP-E)*X
            Y(1:NIN) = DCMPLX(RES(1:NIN)-E0,APP(1:NIN)) * X (1:NIN)

      case (TAAPR)   !  (H+iAPR-E)*X
            Y(1:NIN) = DCMPLX(RES(1:NIN)-E0,APR(1:NIN)) * X (1:NIN)

      case (TMAP)    !  (H-iAP-E)*X
            Y(1:NIN) = DCMPLX(RES(1:NIN)-E0,-AP(1:NIN)) * X (1:NIN)

      case (TMAPP)   !  (H-iAPP-E)*x
            Y(1:NIN) = DCMPLX(RES(1:NIN)-E0,-APP(1:NIN)) * X (1:NIN)
 
      case (TMAPR)   !  (H-iAPR-E)*X
            Y(1:NIN) = DCMPLX(RES(1:NIN)-E0,-APP(1:NIN)) * X (1:NIN)
 
      case default   !default, (H-E)*X
            Y(1:NIN) = (RES(1:NIN)-E0) * X (1:NIN)
      end select


      do level = 1, sF	          
         if (sNDVR .and. (level==sF)) then       
            call H3X_CX(myDim(level),sN(level),OUTHCX, X, tmp)    
         else
      	    if (sDEP(level)) then
                call H1X_DEP_CX (myDim(level), sN(level), myBLK(level),     &
                     H0CX(myH0%mStart(level)),DEPCX(myDEP%mStart(level)), X, tmp)
            else
                call H1X_XYZ_CX (myDim(level), sN(level), myBLK(level),     &
                        H0CX(myH0%mStart(level)), X, tmp)
            end if
         end if

         Y(1:NIN) = Y(1:NIN) + tmp(1:NIN)     
       end do
  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



