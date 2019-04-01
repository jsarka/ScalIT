
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                  Do the real work                       c
!c       Y=Vf*V(f-1)*...*V1*P*V1^T*V2^T*...*Vf^T*X         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine MYPX(nOSB, nPC, N, X, Y)
      integer, intent(IN) :: nOSB, nPC, N
      double precision, intent(IN)  :: X(N)
      double precision, intent(OUT) :: Y(N)

!***********************************************
      integer :: level, i, info, N1
      double precision :: RDE, ADE, DE1, DE2, E0
      double precision :: alpha, beta
      double precision :: Y0(pMax)
!**************************************

      ADE = DABS(sOSBW%mDE);   RDE  = 1.0D0 / ADE
      DE1 =  - ADE;            DE2 =  + ADE  
      beta = 1.0D0 + sOSBW%mbeta       
      alpha = -sOSBW%mbeta*ADE       

      Y0(1:N) = X(1:N)    
      do level=sF, 2, -1
         if (myNode%nodNum(level)>1) then
	    call MG1VTX(level,level-1,nout(level),VOSB(myVOSB%pStart(level)), &
                        nin(level),plen(level-1),Y0)
         else
  	    call MS1VTX(level,level-1,nout(level),VOSB(myVOSB%pStart(level)), &
                        nin(level),blk(level),plen(level-1),Y0)
         end if
      end do

      level = 1
      call VTX_Seq(nin(level), nout(level), blk(level),VOSB(myVOSB%pStart(level)),Y0)

      ! if (id==0) print *, 'nPC=',nPC,' nOSB=',nOSB
      ! if(id==0) print *, 'Node=',ID,',Y0=',Y0(1:plen(1))
      ! if(id==0) print *, 'Node=',ID,',Eig=',EIG0(1:plen(1))
      select case (nPC)
      case (TA0, TAAP0, TMAP0)
          E0 = 0.0D0

      case (TA1)
          E0 = sOSBW%mE0

      case default
          E0 = sOSBW%mE0
      end select

      n1 = plen(1)
      PX_KERNEL: select case (nOSB)
      case (TOSB)     ! 1/Eig,  (Eig-E0)^-1
            Y0(1:N1) = Y0(1:N1)/( EIG0(1:N1) - E0) 

      case (TOSBD1)         ! simple version of OSBD
           do I = 1, N1
              ADE = EIG0(i) - E0
              if (( ADE <= 0.0D0) .and. (ADE >= DE1)) then
                  Y0(i) = -Y0(i)*RDE
              else          
                 if (( ADE <= DE2) .and. (ADE >= 0.0D0)) then
                     Y0(i) = Y0(i)*RDE
                 else             
                     Y0(i) = Y0(i)/ADE
                 end if
              end if
           end do

      case (TOSBD2)      ! (alpha+(1+beta)*DE)^-1
           do I = 1, N1
              ADE = EIG0(i) - E0
              if (( ADE <= 0.0D0) .and. (ADE >= DE1)) then
                   Y0(i) = Y0(i)/((1.0D0+beta)*ADE - alpha)
              else          
                 if (( ADE <= DE2) .and. (ADE >= 0.0D0)) then
                     Y0(i) = Y0(i)/((1.0D0+beta)*ADE + alpha)
                 else             
                     Y0(i) = Y0(i)/ADE
                 end if
              end if
           end do

      case (TOSBW)            
              ! if (hwLen > 0 ) then
            call getWinVector(Y0,WY)       ! WY(1:hwLen) = Y(WIND(1:hwLen))
           
              !  if (id==0) print *, 'initial Vector:', WY
 
            HWTMP(1:hwLen,1:hwLen) = HWR(1:hwLen,1:hwLen)
         
              ! call LAPACK subroutine to get WY=A^-1 * WY
            call DSYTRS('U',hwLen,1,HWTMP,hwLen,WIPIV,WY,hwLen,INFO) 
              !end if

              !if (id==0) print *, 'final Vector A^-1*WY:', WY
          Y0(1:N1) = Y0(1:N1)/( EIG0(1:N1) - E0)

              ! if (id==2) print *, 'N1=',N1,' pind:', pind,  &
              ! ' scnt:',scnt, 'gcnt=',gcnt, 'phwlen=',phwlen
          !if ( hwLen > 0 )   then
  	     do i = 1, phwLen
 	         Y0(pInd(i)) = WY(sCnt(id+1)+i)
             end do
          ! end if

      case default     ! 1/Eig,  (Eig-E0)^-1
            Y0(1:N1) = Y0(1:N1)/( EIG0(1:N1) - E0)

      end select PX_KERNEL
      
      ! if(id==0) print *, 'After kernel, Node=',ID,',Y0=',Y0(1:N1)
      do level=1,sF-1
         if (myNode%nodNUm(level)>1) then
	    call MG1VX(level,level+1,nout(level),VOSB(myVOSB%pStart(level)), &
                        nin(level),plen(level+1),Y0)
         else
  	    call MS1VX(level,level+1,nout(level),VOSB(myVOSB%pStart(level)), &
                        nin(level),blk(level),plen(level+1),Y0)
         end if
      end do

      level = sF
      call VX_Seq(nin(level), nout(level), blk(level),VOSB(myVOSB%pStart(level)),Y0)

      !if(id==0) print *, 'After VX, Node=',ID,',Y0=',Y0(1:N)

      ! Y=Vf*Y0
      Y(1:N)=Y0(1:N)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
