!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            Do the real work                        c
!c               Y = (H-E)*X or H*X                   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 subroutine  MYHX(nType, N, X, Y)    
      integer, intent(IN)    :: nType, N
      double precision,intent(IN)  :: X(N)
      double precision,intent(OUT) :: Y(N)

      double precision :: tmp(pMax), X0(pMax), Y0(pMax)
      integer :: level

      !cccccccccccccccccccccccccccccccccccccccccccccccccccc

      X0(1:N)=X(1:N); Y0(1:N)=0.0D0      

      do level = sF, 2, -1     ! Y0 = Y0 + H*X0, then update X0 and Y0 

         if (sNDVR .and. (level==sF)) then  
            call MG3X_Out(level,level-1,nout(sF),outH,nin(sF),X0,plen(sF-1),Y0)
         else
  	    if (sDEP(level)) then
	       if (myNode%nodNum(level)>1) then
                 call MG2X_Dep(level,level-1,nout(level),H0(myH0%mStart(level))&
                                ,nin(level),DEP(myDEP%pStart(level)),X0,       &
                                plen(level-1),Y0)
               else
                 call MS2X_Dep(level,level-1,nout(level),H0(myH0%mStart(level))&
                               ,nin(level),blk(level),DEP(myDEP%pStart(level)),&
                                X0, plen(level-1),Y0)
               end if
            else             
  	       if (myNode%nodNum(level)>1) then
                 call MG1X_XYZ(level,level-1,nout(level),H0(myH0%mStart(level))&
                                ,nin(level),X0,plen(level-1),Y0 )
               else
                 call MS1X_XYZ(level,level-1,nout(level),H0(myH0%mStart(level))&
                                ,nin(level),blk(level),X0,plen(level-1),Y0 ) 
               end if
            end if
         end if

      end do
      
      !  For level = 1, no coordinate dependence,
      !  we assume that #of nodes < # of blocks.
      level = 1
      call H1X_XYZ_Seq(nin(level),nout(level),blk(level),H0(myH0%mStart(level)),  &
                   X0,tmp)

      select case(nType)
          case (TA0) 
              Y0(1:plen(1)) = tmp(1:plen(1))+ResSeq(1:plen(1))*X0(1:plen(1))  &
                              + Y0(1:plen(1))
          case (TA1)
              Y0(1:plen(1)) = tmp(1:plen(1))+( ResSeq(1:plen(1))-sOSBW%mE0)   &
                                * X0(1:plen(1)) + Y0(1:plen(1))
          case default
              Y0(1:plen(1)) = tmp(1:plen(1))+( ResSeq(1:plen(1))-sOSBW%mE0)   &
                                * X0(1:plen(1)) + Y0(1:plen(1))
      end select        
      
      call updateX(1,sF,Y0,Y)

   end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
