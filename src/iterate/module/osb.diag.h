!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function BlockDiag()

    BlockDiag = myDiag(myDiagInit, myDiagInitCX)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function BlockDiagUser (userFunc, userFuncCX)
    external :: userFunc, userFuncCX

    BlockDiagUser = myDiag(userFunc, userFuncCX)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine myDiagInit(level)
     integer, intent(IN) :: level
     
     if (sNDVR .AND. (level==sF)) then
        HOSB(myHOSB%mStart(sF):myHOSB%mEnd(sF))=OUTH(1:sN(sF)*myLen)  
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
        HOSBCX(myHOSB%mStart(sF):myHOSB%mEnd(sF))=OUTHCX(1:sN(sF)*myLen)    
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
logical function myDiag(userFunc, userFuncCX)
   external :: userFunc, userFuncCX

   integer  :: level, num, info, i, sumOff
   double precision :: db
   double complex ::  dbcx
   logical :: saveHB   

   myDiag = .false.

   saveHB = ((sHOSB>0) .AND. (.NOT. sST))

   if (sCX) then
      inquire(IOLENGTH=num) dbcx
   else
      inquire(IOLENGTH=num) db
   end if  

   if (saveHB) then
      open(99, file=HOSBFile, status='Replace', FORM='UNFORMATTED', &
           access='direct', recl=num, IOSTAT=info)
      if (info/=0) return
   end if
   
   if (sCX) then

      do level = sF, 1, -1

          print *
          write(*, 10) level,sN(level),myBlk(level)
          write(*, 11) myDim(level),sDep(level)

          if (.NOT. sST) then
             if ( .NOT. allocHOSBCX(sN(level)*myLen))  then
	        print *, ' Error in allocating HOSB at level ', level
                exit 
             end if
          end if

          call userFuncCX(level)
          call HosbDiagDX(sBJ%mMax, sBJ%mTol, myBlk(level), sN(level), &
	                  myDim(level), HOSBCX(myHOSB%mStart(level)),  &
                          VOSB(myVOSB%mStart(level)), EIG0)

         if (saveHB) then
            sumOff = (Sum(sN(1:level))-sN(level))*myLen
  	    do i = 1, sN(level)*myLen
 	        num = sumOff + i
 	        write(99,REC=num) HOSBCX(i)
            end do
             
         end if
      end do

      if ((.NOT. sST).AND.(allocated(HOSBCX)))  deallocate(HOSBCX)

   else

      do level = sF, 1, -1

          print *
          write(*, 10) level,sN(level),myBlk(level)
          write(*, 11) myDim(level),sDep(level)

          if (.NOT. sST) then
             if ( .NOT. allocHOSB(sN(level)*myLen))  then
	        print *, ' Error in allocating HOSB at level ', level
  	        exit 
             end if
          end if

          call userFunc(level)

          call HosbDiag(sBJ%mMax, sBJ%mTol, myBlk(level), sN(level), &
	                myDim(level), HOSB(myHOSB%mStart(level)),    &
                        VOSB(myVOSB%mStart(level)), EIG0)

          if (saveHB) then
             sumOff = (Sum(sN(1:level))-sN(level))*myLen
    	     do i = 1, sN(level)*myLen
 	         num = sumOff + i
 	         write(99,REC=num) HOSB(i)
             end do
          end if
      end do

      if ((.NOT. sST).AND.(allocated(HOSB)))  deallocate(HOSB)

   end if

   if (saveHB) close(99)
   
   myDiag = .true.

 10  FORMAT('  BJ Diag.: Level:',I4,'.   Grid Size:',I6, &
            '   # of Blocks:', I6)
 11  FORMAT('  Dimension Size:',I10,'.   Coord. Dep:',L4)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
