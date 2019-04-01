
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                  Do the real work                       c
!c       Y=Vf*V(f-1)*...*V1*P*V1^T*V2^T*...*Vf^T*X         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine MYPX_DX(nOSB, nPC, N, X, Y)
      integer, intent(IN) :: nOSB,nPC,N
      double complex, intent(IN)  :: X(N)
      double complex, intent(OUT) :: Y(N)

!***********************************************
      integer :: level, i, info, osbwcnt, N1
      double precision :: RDE, ADE, DE1, DE2, E0, DE
      double precision :: alpha, beta
double complex :: X0(pMax), Y0(pMax), tmp
!**************************************

      DE = DABS(sOSBW%mDE);   RDE  = 1.0D0 / DE
      DE1 =  - DE;            DE2 =  + DE  

      Y0(1:N) = X(1:N)    
      do level=sF, 2, -1
         if (myNode%nodNum(level)>1) then
	    call MG1VTX_DX(level,level-1,nout(level),VOSB(myVOSB%pStart(level)), &
                        nin(level),plen(level-1),Y0)
         else
  	    call MS1VTX_DX(level,level-1,nout(level),VOSB(myVOSB%pStart(level)), &
                        nin(level),blk(level),plen(level-1),Y0)
         end if
      end do

      level = 1
      call VTX_DX_Seq(nin(level), nout(level), blk(level),VOSB(myVOSB%pStart(level)),Y0)

      select case (nPC)
      case (TA0,TAAP0,TMAP0)
          E0 = 0.0D0

      case (TA1)
          E0 = sOSBW%mE0

      case default
          E0 = sOSBW%mE0
      end select

      N1 = plen(1)
      PX_KERNEL: select case (nOSB)
      case (TOSB)     ! 1/Eig,  (Eig-E0)^-1
         select case (nPC)
         case (TA0)
            Y0(1:N1) = Y0(1:N1)/EIG0(1:N1) 
         case (TA1)
  	    Y0(1:N1) = Y0(1:N1)/(EIG0(1:N1)-E0) 
         case (TAAP0)
 	    Y0(1:N1) = Y0(1:N1)/DCMPLX(EIG0(1:N1),AP(1:N1)) 
         case (TMAP0)
 	    Y0(1:N1) = Y0(1:N1)/DCMPLX(EIG0(1:N1),-AP(1:N1)) 
         case (TAAP)
 	    Y0(1:N1) = Y0(1:N1)/DCMPLX(EIG0(1:N1)-E0,AP(1:N1)) 
         case (TAAPP) 
            Y0(1:N1) = Y0(1:N1)/DCMPLX(EIG0(1:N1)-E0,APP(1:N1))
         case (TAAPR) 
            Y0(1:N1) = Y0(1:N1)/DCMPLX(EIG0(1:N1)-E0,APR(1:N1))
         case (TMAP) 
            Y0(1:N1) = Y0(1:N1)/DCMPLX(EIG0(1:N1)-E0,-AP(1:N1))
         case (TMAPP) 
            Y0(1:N1) = Y0(1:N1)/DCMPLX(EIG0(1:N1)-E0,-APP(1:N1))
         case (TMAPR) 
            Y0(1:N1) = Y0(1:N1)/DCMPLX(EIG0(1:N1)-E0,-APR(1:N1))
         case default
            Y0(1:N1) = Y0(1:N1)/(EIG0(1:N1)-E0)
         end select

      case (TOSBD1)         ! simple version of OSBD
           do I = 1, N1
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
             Y0(i) = Y0(i)/tmp
           end do

      case (TOSBD2)      ! (alpha+(1+beta)*DE)^-1
         beta  = sOSBW%mBeta+1.0D0;
         alpha = -sOSBW%mBeta*DE

         do I = 1, N1
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
             Y0(i) = Y0(i)/tmp
         end do 

!           do I = 1, N1
!              ADE = EIG0(i) - E0
!              if (( ADE <= 0.0D0) .and. (ADE >= DE1)) then
!                   Y0(i) = Y0(i)/((1.0D0+beta)*ADE - alpha)
!              else          
!                 if (( ADE <= DE2) .and. (ADE >= 0.0D0)) then
!                     Y0(i) = Y0(i)/((1.0D0+beta)*ADE + alpha)
!                 else             
!                     Y0(i) = Y0(i)/ADE
!                 end if
!              end if
!           end do
!
!      case (TOSBW)            
!         if (hwLen > 0 ) then
!            call getWinVector_CX(Y0,WYCX)       ! WY(1:hwLen) = Y(WIND(1:hwLen))
!            
!            HWTMPCX(1:hwLen,1:hwLen) = HWRCX(1:hwLen,1:hwLen)
!         
!            ! call LAPACK subroutine to get WY=A^-1 * WY
!  	    if ((nPC==TA0) .OR. (nPC==TA1)) then
!                call ZHETRS('U',hwLen,1,HWTMPCX,hwLen,WIPIV,WYCX,hwLen,INFO) 
!            else
!                call ZGETRS('N',hwLen,1,HWTMPCX,hwLen,WIPIV,WYCX,hwLen,INFO) 
!            end if
!
!         end if
!
!          Y0(1:N1) = Y0(1:N1)/( EIG0(1:N1) - E0)
!
!          if ( hwLen > 0 )   then
!  	     do i = 1, phwLen
! 	         Y0(pInd(i)) = WYCX(sCnt(id+1)+i)
!             end do
!           end if
!
       case (TOSBW)         
         call getWinVector_CX(Y0,WYCX)       ! WY(1:hwLen) = Y(WIND(1:hwLen))
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
             Y0(1:N1) = Y(1:N1)/EIG0(1:N1)
         case (TA1)
             Y0(1:N1) = Y(1:N1)/(EIG0(1:N1)-E0)
         case (TAAP,TAAP0)
             Y0(1:N1) = Y(1:N1)/DCMPLX(EIG0(1:N1)-E0, AP(1:N1))
         case (TAAPP)
             Y0(1:N1) = Y(1:N1)/DCMPLX(EIG0(1:N1)-E0, APP(1:N1))
         case (TAAPR)
             Y0(1:N1) = Y(1:N1)/DCMPLX(EIG0(1:N1)-E0, APR(1:N1))
         case (TMAP,TMAP0)
             Y0(1:N1) = Y(1:N1)/DCMPLX(EIG0(1:N1)-E0, -AP(1:N1))
         case (TMAPP)
             Y0(1:N1) = Y(1:N1)/DCMPLX(EIG0(1:N1)-E0, -APP(1:N1))
         case (TMAPR)
             Y0(1:N1) = Y(1:N1)/DCMPLX(EIG0(1:N1)-E0, -APR(1:N1))
         case default
             Y0(1:N1) = Y(1:N1)/(EIG0(1:N1)-E0)
         end select
 
         do i = 1, phwLen             ! Restore Y(index(I))
	    Y0(pInd(i)) = WYCX(sCnt(id+1)+i)
         end do          

      case default
         select case (nPC)
         case (TA0)
             Y(1:N1) = Y(1:N1)/EIG0(1:N1)
         case (TA1)
             Y(1:N1) = Y(1:N1)/(EIG0(1:N1)-E0)
         case (TAAP,TAAP0)
             Y(1:N1) = Y(1:N1)/DCMPLX(EIG0(1:N1)-E0, AP(1:N1))
         case (TAAPP)
             Y(1:N1) = Y(1:N1)/DCMPLX(EIG0(1:N1)-E0, APP(1:N1))
         case (TAAPR)
             Y(1:N1) = Y(1:N1)/DCMPLX(EIG0(1:N1)-E0, APR(1:N1))
         case (TMAP,TMAP0)
             Y(1:N1) = Y(1:N1)/DCMPLX(EIG0(1:N1)-E0, -AP(1:N1))
         case (TMAPP)
             Y(1:N1) = Y(1:N1)/DCMPLX(EIG0(1:N1)-E0, -APP(1:N1))
         case (TMAPR)
             Y(1:N1) = Y(1:N1)/DCMPLX(EIG0(1:N1)-E0, -APR(1:N1))
         case default
             Y(1:N1) = Y(1:N1)/(EIG0(1:N1)-E0)
         end select        
      end select PX_KERNEL


      do level=1,sF-1
         if (myNode%nodNUm(level)>1) then
	    call MG1VX_DX(level,level+1,nout(level),VOSB(myVOSB%pStart(level)), &
                        nin(level),plen(level+1),Y0)
         else
  	    call MS1VX_DX(level,level+1,nout(level),VOSB(myVOSB%pStart(level)), &
                        nin(level),blk(level),plen(level+1),Y0)
         end if
      end do

      level = sF
      call VX_DX_Seq(nin(level), nout(level), blk(level),VOSB(myVOSB%pStart(level)),Y0)

      ! Y=Vf*Y0
      Y(1:N)=Y0(1:N)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
