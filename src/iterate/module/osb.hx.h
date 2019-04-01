!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  HX part of OSB     :                             c
!c  Subroutines:                                     c
!c  HX    (NIN, X, Y):  Y = H * X                    c
!c  EHX   (NIN, X, Y):  Y = ( H - E) * X             c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccc
!c             Y = H * X                    c
!c        X, Y should be different          c
!cccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine  HX(NIN, X, Y)
      integer, intent(IN)    :: NIN
      double precision,intent(IN)  :: X(NIN)
      double precision,intent(OUT) :: Y(NIN)

      call MYHX(sHC, NIN, X, Y)

 end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccc
!c             Y = H * X                    c
!c        X, Y should be different          c
!cccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine  HijX(NIN, X, Y)
      integer, intent(IN)    :: NIN
      double precision,intent(IN)  :: X(NIN)
      double precision,intent(OUT) :: Y(NIN)

      call MYHX(TA0, NIN, X, Y)

 end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine  HX0(NIN, X, Y)
      integer, intent(IN)    :: NIN
      double precision,intent(IN)  :: X(NIN)
      double precision,intent(OUT) :: Y(NIN)

      call MYHX(TA0, NIN, X, Y)

 end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Y = ( H - E ) * X                       c
!c    X and Y should not be the same.                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  EHX(NIN, X, Y)    
      integer, intent(IN)    :: NIN     
      double precision,intent(IN)  :: X(NIN)
      double precision,intent(OUT) :: Y(NIN)

      call MYHX(TA1, NIN, X, Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Do the real work                        c
!c               Y = (H-E)*X or H*X                   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine  MYHX(nType, NIN, X, Y)    
      integer, intent(IN)    :: nType
      integer, intent(IN)    :: NIN     
      double precision,intent(IN)  :: X(NIN)
      double precision,intent(OUT) :: Y(NIN)

      double precision :: tmp(NIN)
      integer :: level

      !cccccccccccccccccccccccccccccccccccccccccccccccccccc
       
      select case(nType)
          case (TA0) 
               Y(1:NIN) = RES(1:NIN) * X (1:NIN)
          case (TA1)
                Y(1:NIN) = ( RES(1:NIN) - sOSBW%mE0 ) * X (1:NIN)
          case default
                Y(1:NIN) = ( RES(1:NIN) - sOSBW%mE0 ) * X (1:NIN)
      end select   


      do level = 1, sF	          
         if (sNDVR .and. (level==sF)) then          
            call H3X(myDim(sF),sN(sF),OUTH, X, tmp)    
	 else
    	    if (sDEP(level)) then
                call H1X_DEP (myDim(level), sN(level), myBlk(level),      &
                          H0(myH0%mStart(level)),DEP(myDEP%mStart(level)), &
                          X, tmp)
            else
                call H1X_XYZ (myDim(level), sN(level), myBLK(level),       &
                        H0(myH0%mStart(level)), X, tmp)
            end if
         end if

         Y(1:NIN) = Y(1:NIN) + tmp(1:NIN)     

      end do

   end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

