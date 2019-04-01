!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Do the real work                        c
!c               Y = (H-E)*X or H*X                   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine  MYHX_CX(nType, N, X, Y)    
      integer, intent(IN)    :: nType, N
      double complex,intent(IN)  :: X(N) 
      double complex,intent(OUT) :: Y(N)

      double complex   :: tmp(pMax), X0(pMax), Y0(pMax)
      integer :: level


      !cccccccccccccccccccccccccccccccccccccccccccccccccccc

      X0(1:N)=X(1:N); Y0(1:N)=0.0D0      

      do level = sF, 2, -1     ! Y0 = Y0 + H*X0, then update X0 and Y0 

         if (sNDVR .and. (level==sF)) then          
            call MG3X_Out_CX(level,level-1,nout(sF),outHCX,nin(sF),X0, &
                             plen(sF-1),Y0)
         else
  	    if (sDEP(level)) then
	       if (myNode%nodNum(level)>1) then
                  call MG2X_Dep_CX(level,level-1,nout(level),          &
                                H0CX(myH0%mStart(level)),nin(level),   &
                                DEPCX(myDEP%pStart(level)),X0,         &
                                plen(level-1),Y0)
               else
                  call MS2X_Dep_CX(level,level-1,nout(level),          &
                                H0CX(myH0%mStart(level)),nin(level),   &
                                blk(level),DEPCX(myDEP%pStart(level)), &
                                X0, plen(level-1),Y0)
               end if
            else             
  	       if (myNode%nodNum(level)>1) then
                  call MG1X_XYZ_CX(level,level-1,nout(level),          &
                                H0CX(myH0%mStart(level)),nin(level),   &
                                X0,plen(level-1),Y0 )
               else
                  call MS1X_XYZ_CX(level,level-1,nout(level),          &
                                H0CX(myH0%mStart(level)),nin(level),   &
                                blk(level),X0,plen(level-1),Y0 ) 
               end if
            end if
         end if
      end do
      
      ! level=1, no coordinate dependence,and #of nodes < #of blocks.
      level = 1
      call H1X_XYZ_CX_Seq(nin(level),nout(level),blk(level),               &
                      H0CX(myH0%mStart(level)),X0,tmp)

      select case(nType)
          case (TA0)     !  H*X
               Y0(1:plen(1)) = tmp(1:plen(1))+ResSeq(1:plen(1))*X0(1:plen(1))  &
                               + Y0(1:plen(1))
          case (TA1)     !  (H-E)*X
               Y0(1:plen(1)) = tmp(1:plen(1))+(ResSeq(1:plen(1))-sOSBW%mE0)    &
                                *X0(1:plen(1)) + Y0(1:plen(1))

          case (TAAP0)   !  (H+iAP-E)*X        
               Y0(1:plen(1))=tmp(1:plen(1))+DCMPLX(ResSeq(1:plen(1)),          &
                        AP(1:plen(1)))* X0(1:plen(1)) + Y0(1:plen(1))

          case (TMAP0)   !  (H-iAP)*X        
               Y0(1:plen(1))=tmp(1:plen(1))+DCMPLX(ResSeq(1:plen(1)),          &
                        -AP(1:plen(1)))* X0(1:plen(1)) + Y0(1:plen(1))

          case (TAAP)    !  (H+iAP-E)*X        
               Y0(1:plen(1))=tmp(1:plen(1))+DCMPLX(ResSeq(1:plen(1))-sOSBW%mE0,&
                        AP(1:plen(1)))* X0(1:plen(1)) + Y0(1:plen(1))

          case (TAAPP)   !  (H+iAPP-E)*X
               Y0(1:plen(1))=tmp(1:plen(1))+DCMPLX(ResSeq(1:plen(1))-sOSBW%mE0,&
                        APP(1:plen(1)))* X0(1:plen(1)) + Y0(1:plen(1))

          case (TAAPR)   !  (H+iAPR-E)*X
               Y0(1:plen(1))=tmp(1:plen(1))+DCMPLX(ResSeq(1:plen(1))-sOSBW%mE0,&
                        APR(1:plen(1)))* X0(1:plen(1)) + Y0(1:plen(1))

          case (TMAP)    !  (H-iAP-E)*X
               Y0(1:plen(1))=tmp(1:plen(1))+DCMPLX(ResSeq(1:plen(1))-sOSBW%mE0,&
                        -AP(1:plen(1)))* X0(1:plen(1)) + Y0(1:plen(1))

          case (TMAPP)   !  (H-iAPP-E)*x
               Y0(1:plen(1))=tmp(1:plen(1))+DCMPLX(ResSeq(1:plen(1))-sOSBW%mE0,&
                        -APP(1:plen(1)))* X0(1:plen(1)) + Y0(1:plen(1))
 
          case (TMAPR)   !  (H-iAPR-E)*X
               Y0(1:plen(1))=tmp(1:plen(1))+DCMPLX(ResSeq(1:plen(1))-sOSBW%mE0,&
                        -APR(1:plen(1)))* X0(1:plen(1)) + Y0(1:plen(1))

          case default
               Y0(1:plen(1)) = tmp(1:plen(1))+( ResSeq(1:plen(1))-sOSBW%mE0)   &
                                *X0(1:plen(1)) + Y0(1:plen(1))
      end select        

      call updateX_CX(1,sF,Y0,Y)

  end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
