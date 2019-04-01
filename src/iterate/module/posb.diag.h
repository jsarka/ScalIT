!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function BlockDiag()

    BlockDiag = myDiag(myDiagInit, myDiagInitCX)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function BlockDiagUser(userFunc, userFuncCX)
    external :: userFunc, userFuncCX

    BlockDiagUser = myDiag(userFunc, userFuncCX)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine myDiagInit(level)
     integer, intent(IN) :: level

     if (sNDVR .AND. (level==sF)) then
        HOSB(myHOSB%pStart(sF):myHOSB%pEnd(sF))=OUTH(1:myHOSB%pSize(sF))
     else
        if (sDep(level)) then
            call myHOSBDepInit(level)
        else
   	    call myHOSBInit(level)
        end if
     end if

end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine myDiagInitCX(level)
     integer, intent(IN) :: level

     if (sNDVR .AND. (level==sF)) then
        HOSBCX(myHOSB%pStart(sF):myHOSB%pEnd(sF))=OUTHCX(1:myHOSB%pSize(sF))    
     else
        if (sDep(level)) then
           call myHOSBDepInit_CX(level)
        else
   	   call myHOSBInit_CX(level)
        end if
     end if

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function myDiag(initFunc, initFuncCX)
   external :: initFunc, initFuncCX

   integer  :: level,  num, ierr, info, fh
   logical  :: sSH
   double precision, allocatable  :: APE0(:)
   double precision :: ct1, ct2, ct3, ct4, MPI_WTime

   myDiag = .false.

   allocate(APE0(pMax),stat=info)  
   if (info/=0) then
        write(*,5) ' E0', level, id;     return
   end if

   sSH = ((.NOT.sST).AND.(sHOSB>0))

   EIG0(1:myRES%pSize(sF)) = RES(1:myRES%pSize(sF))

   if (sSH)   call MPI_FILE_Open(MPI_COMM_WORLD,fHOSB,  &
            MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)

   if (sCX) then   
  
     do level = sF, 1, -1

          if (id==rootID) then
             print *
             write(*, 10) level,sN(level),blk(level)
             write(*, 11) nin(level),sDep(level)
          end if

          if (.NOT. sST) then
             if ( .NOT. allocHOSBCX(myHOSB%pSize(level)))  then
  	        write(*,5) ' HOSBCx', level, id
                deallocate(APE0);     exit  
             end if
          end if

          call initFuncCX(level)

          if (myNode%nodNum(level)>1) then
              call HosbDiag_DX_MPI(myNode%commID(level),sBJ%mMax, sBJ%mTol, &
                          blk(level), nout(level), nin(level),              &
                          HOSBCX(myHOSB%pStart(level)),                     &
                          VOSB(myVOSB%pStart(level)), APE0, ierr)
          else
              call HosbDiagDX_Seq(sBJ%mMax,sBJ%mTol,blk(level),nout(level),   &
	                  nin(level), HOSBCX(myHOSB%pStart(level)),       &
                          VOSB(myVOSB%pStart(level)), APE0)
          end if

          ! update RES
          if (level>1) then
              call UpdateX(level, level-1, APE0, EIG0)
          else
  	      Eig0(1:myRES%pSize(1)) = APE0(1:myRES%pSize(1))
          end if

          if (sSH) then    ! store data
	     if (myNode%nodNum(level)>1) then    ! store grid data
 	        call MWriteDataGrid_CX(fh,cxSize,myHOSB%gPos(level),     &
                      myconf%gDim(level),sN(level),myData%pDim(level),   &
                           HOSBCX(myHOSB%pStart(level)), ierr)
             else           ! store sequential data
	        call MWriteData_CX(fh,cxSize,myHOSB%gPos(level),         &
		      myHOSB%pSize(level),HOSBCX(myHOSB%pStart(level)),ierr)
             end if
          end if
      end do

      if ((.NOT. sST).AND.(allocated(HOSBCX)))  deallocate(HOSBCX)

   else
      do level = sF, 1, -1
          if (id==rootID) then
             print *
             write(*, 10) level,sN(level),blk(level)
             write(*, 11) nin(level),sDep(level)
             ct1 = MPI_WTIME()
          end if

          if (.NOT. sST) then
             if ( .NOT. allocHOSB(myHOSB%pSize(level)))  then
 	         write(*,5) ' HOSB', level, id
                 deallocate(APE0);        exit 
             end if
          end if

          call initFunc(level)

          if (myNode%nodNum(level)>1) then          
             call HosbDiag_MPI(myNode%commID(level),sBJ%mMax,sBJ%mTol, &
                             blk(level), nout(level), nin(level),      &
                             HOSB(myHOSB%pStart(level)),               &
                             VOSB(myVOSB%pStart(level)), APE0, ierr, myNode%id)
          else
             call HosbDiag_Seq(sBJ%mMax,sBJ%mTol,blk(level),nout(level),   &
	                nin(level), HOSB(myHOSB%pStart(level)),        &
                        VOSB(myVOSB%pStart(level)), APE0, myNode%id)
          end if   

          if (id==rootID) then
              ct2 = MPI_WTIME()
              print *
              print *, '   MPI WTime for HOSB compute (sec):',ct2-ct1
          end if

          if (level>1)  then
              call UpdateX(level, level-1, APE0, EIG0)
          else
  	      Eig0(1:myRES%pSize(1)) = APE0(1:myRES%pSize(1))
          end if

          if (sSH) then  ! store HOSB data
              if (myNode%nodNum(level)>1) then    ! store grid data
                 call MWriteDataGrid(fh,dbSize,myHOSB%gPos(level),         &
                   myconf%gDim(level),sN(level),HOSB(myHOSB%pStart(level)),ierr)
              else           ! store sequential data
                 call MWriteData(fh,dbSize,myHOSB%gPos(level),             &
           	   myHOSB%pSize(level),HOSB(myHOSB%pStart(level)),ierr)
              end if
          end if          

          if (id==rootID) then
              ct3 = MPI_WTIME()
              print *, '   MPI WTime for HOSB update  (sec):',ct3-ct2
              print *
          end if

      end do

      if ((.NOT. sST).AND.(allocated(HOSB)))  deallocate(HOSB)

   end if

   if (sSH)  call MPI_FILE_CLOSE(fh, ierr)

   deallocate(APE0)
   myDiag = .true.

  5  format(' Error in allocating ', A,' at layer = ',I5,' with node id = ',I6)
 10  format('  BJ Diag.: Level:',I4,'.   Grid Size:',I6,'   # of Blocks:',I6)
 11  format('  Dimension Size:',I10,'.   Coord. Dep:',L4)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
