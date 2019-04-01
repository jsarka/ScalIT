!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                  Do the real work                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine MYPX_DX(nOSB, nPC, NIN, X, Y)
      integer, intent(IN) :: nOSB, nPC, NIN
      double complex, intent(IN)  :: X(NIN)
      double complex, intent(OUT) :: Y(NIN)

!***********************************************
      integer :: level, i, info
      double precision :: DE, RDE, ADE, E0
      double precision :: alpha, beta
      double complex :: tmp
!**************************************

      Y(1:NIN) = X(1:NIN)

      do level = sF, 1, -1
          call VTX_DX(myDim(level), sN(level), myBLK(level), &
                      VOSB(myVOSB%mStart(level)),Y)
      end do

      select case (nPC)
      case (TA0,TAAP0,TMAP0)
         E0 = 0.0D0

      case default
         E0 = sOSBW%mE0
      end select

      DE = DABS(sOSBW%mDE);  RDE=1.0D0/DE
      
!***************************************************************
      KERNEL_PXCX : select case (nOSB)
      case (TOSB)
         select case (nPC)
         case (TA0) 
             Y(1:NIN) = Y(1:NIN)/EIG0(1:NIN)
         case (TA1) 
             Y(1:NIN) = Y(1:NIN)/(EIG0(1:NIN)-E0)
         case (TAAP0) 
             Y(1:NIN) = Y(1:NIN)/DCMPLX(EIG0(1:NIN), AP(1:NIN))
         case (TMAP0) 
             Y(1:NIN) = Y(1:NIN)/DCMPLX(EIG0(1:NIN), -AP(1:NIN))
         case (TAAP) 
             Y(1:NIN) = Y(1:NIN)/DCMPLX(EIG0(1:NIN)-E0, AP(1:NIN))
         case (TAAPP) 
             Y(1:NIN) = Y(1:NIN)/DCMPLX(EIG0(1:NIN)-E0, APP(1:NIN))
         case (TAAPR) 
             Y(1:NIN) = Y(1:NIN)/DCMPLX(EIG0(1:NIN)-E0, APR(1:NIN))
         case (TMAP) 
             Y(1:NIN) = Y(1:NIN)/DCMPLX(EIG0(1:NIN)-E0, -AP(1:NIN))
         case (TMAPP) 
             Y(1:NIN) = Y(1:NIN)/DCMPLX(EIG0(1:NIN)-E0, -APP(1:NIN))
         case (TMAPR) 
             Y(1:NIN) = Y(1:NIN)/DCMPLX(EIG0(1:NIN)-E0, -APR(1:NIN))
         case default
             Y(1:NIN) = Y(1:NIN)/(EIG0(1:NIN)-E0)
         end select

      case (TOSBD1)         
         do I = 1, NIN
             ADE = EIG0(i)-E0
             if ( ABS(ADE) <= DE ) then
	        if (ADE<0.0D0) then
                   ADE = -DE
                else
                   ADE = DE
                end if
             end if
             select case(nPC)
             case (TA0, TA1)
                 tmp = ADE             
             case (TAAP, TAAP0) 
	         tmp = DCMPLX(ADE, AP(i))
             case (TAAPP) 
	         tmp = DCMPLX(ADE, APP(i))
             case (TAAPR) 
	         tmp = DCMPLX(ADE, APR(i))
             case (TMAP, TMAP0) 
	         tmp = DCMPLX(ADE, -AP(i))
             case (TMAPP) 
	         tmp = DCMPLX(ADE, -APP(i))
             case (TMAPR) 
	         tmp = DCMPLX(ADE, -APR(i))
             case default
                 tmp = ADE 
             end select
             Y(i) = Y(i)/tmp
         end do             

      case (TOSBD2)
         beta  = sOSBW%mBeta+1.0D0;
         alpha = -sOSBW%mBeta*DE

         do I = 1, NIN
             ADE = EIG0(i)-E0
             if ( ABS(ADE) <= DE ) then
	        if (ADE<0) then
		    ADE = alpha + beta*ADE
                else
		    ADE = alpha + beta*ADE
                end if 	      
             end if             
             select case(nPC)                
             case (TA0, TA1)
                 tmp = ADE             
             case (TAAP,TAAP0) 
	         tmp = DCMPLX(ADE, AP(i))
             case (TAAPP) 
	         tmp = DCMPLX(ADE, APP(i))
             case (TAAPR) 
	         tmp = DCMPLX(ADE, APR(i))
             case (TMAP,TMAP0) 
	         tmp = DCMPLX(ADE, -AP(i))
             case (TMAPP) 
	         tmp = DCMPLX(ADE, -APP(i))
             case (TMAPR) 
	         tmp = DCMPLX(ADE, -APR(i))
             case default
                 tmp = ADE 
             end select
             Y(i) = Y(i)/tmp
         end do 

      case (TOSBW)         
         do i = 1, hwLen
            WYCX(i) = Y(WIND(i))
         end do
       
	 HWTMPCX(1:hwLen,1:hwLen) = HWRCX(1:hwLen,1:hwLen)

 	 select case (nPC)
	 case (TA0,TA1)
            call ZHETRS('U', hwLen, 1, HWTMPCX, hwLen, WIPIV,       &
                       WYCX, hwLen, INFO)

	 case (TAAP, TAAPP, TAAPR, TMAP, TMAPP, TMAPR, TAAP0, TMAP0)
 
            call ZGETRS('N', hwLen, 1, HWTMPCX, hwLen, WIPIV,       &
                         WYCX, hwLen, INFO)
	 case default
            call ZHETRS('U', hwLen, 1, HWTMPCX, hwLen, WIPIV,       &
                       WYCX, hwLen, INFO)
         end select


         select case (nPC)
         case (TA0)
             Y(1:NIN) = Y(1:NIN)/EIG0(1:NIN)
         case (TA1)
             Y(1:NIN) = Y(1:NIN)/(EIG0(1:NIN)-E0)
         case (TAAP,TAAP0)
             Y(1:NIN) = Y(1:NIN)/DCMPLX(EIG0(1:NIN)-E0, AP(1:NIN))
         case (TAAPP)
             Y(1:NIN) = Y(1:NIN)/DCMPLX(EIG0(1:NIN)-E0, APP(1:NIN))
         case (TAAPR)
             Y(1:NIN) = Y(1:NIN)/DCMPLX(EIG0(1:NIN)-E0, APR(1:NIN))
         case (TMAP,TMAP0)
             Y(1:NIN) = Y(1:NIN)/DCMPLX(EIG0(1:NIN)-E0, -AP(1:NIN))
         case (TMAPP)
             Y(1:NIN) = Y(1:NIN)/DCMPLX(EIG0(1:NIN)-E0, -APP(1:NIN))
         case (TMAPR)
             Y(1:NIN) = Y(1:NIN)/DCMPLX(EIG0(1:NIN)-E0, -APR(1:NIN))
         case default
             Y(1:NIN) = Y(1:NIN)/(EIG0(1:NIN)-E0)
         end select
 
         do i = 1, hwLen             ! Restore Y(index(I))
	    Y(WIND(i)) = WYCX(i)
         end do          

      case default
         select case (nPC)
         case (TA0)
             Y(1:NIN) = Y(1:NIN)/EIG0(1:NIN)
         case (TA1)
             Y(1:NIN) = Y(1:NIN)/(EIG0(1:NIN)-E0)
         case (TAAP,TAAP0)
             Y(1:NIN) = Y(1:NIN)/DCMPLX(EIG0(1:NIN)-E0, AP(1:NIN))
         case (TAAPP)
             Y(1:NIN) = Y(1:NIN)/DCMPLX(EIG0(1:NIN)-E0, APP(1:NIN))
         case (TAAPR)
             Y(1:NIN) = Y(1:NIN)/DCMPLX(EIG0(1:NIN)-E0, APR(1:NIN))
         case (TMAP,TMAP0)
             Y(1:NIN) = Y(1:NIN)/DCMPLX(EIG0(1:NIN)-E0, -AP(1:NIN))
         case (TMAPP)
             Y(1:NIN) = Y(1:NIN)/DCMPLX(EIG0(1:NIN)-E0, -APP(1:NIN))
         case (TMAPR)
             Y(1:NIN) = Y(1:NIN)/DCMPLX(EIG0(1:NIN)-E0, -APR(1:NIN))
         case default
             Y(1:NIN) = Y(1:NIN)/(EIG0(1:NIN)-E0)
         end select
         
      end select KERNEL_PXCX
!*************************************************************************

      do level = 1, sF
         call VX_DX(myDim(level), sN(level), myBLK(level),    &
                    VOSB(myVOSB%mStart(level)), Y)
      end do

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

