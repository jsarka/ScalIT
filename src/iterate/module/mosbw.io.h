!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine saveHW()
   integer :: ierr

   if (id == rootID) then
      open(99,FILE=fHW, STATUS='REPLACE',FORM='UNFORMATTED')

       if (sCX) then
	  write(99) HWCX	        
       else
	  write(99) HW                         
       end if             

       close(99)
   end if

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine loadHW()
   integer :: ierr

   if (id == rootID) then
      open(99,FILE=fHW, STATUS='OLD',FORM='UNFORMATTED')

       if (sCX) then
	  read(99) HWCX	        
       else
	  read(99) HW                         
       end if             

       close(99)
   end if

   if (sCX) then
      call MPI_BCAST(HWCX,hwLen**2,MPI_DOUBLE_COMPLEX,rootID,    &
                         MPI_COMM_WORLD,ierr)
   else
      call MPI_BCAST(HW,hwLen**2,MPI_DOUBLE_PRECISION,rootID,    &
                         MPI_COMM_WORLD,ierr)
   end if

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function saveHWLocal()

    integer :: info
 
    saveHWLocal = .false.

    if (srMode) then
          open(99,FILE=fHW, STATUS='REPLACE', FORM='UNFORMATTED', &
               IOSTAT=info)
          if (info /= 0) return
          write(99) sOSBW%mCnt,sOSBW%mE0,sOSBW%mDE,sOSBW%mBeta,hwLen
	  if (sCX) then
	       write(99) HWCX
          else
               write(99) HW
          end if
          close(99)
          saveHWLocal = .true.
      else
          open(99,FILE=fHW, STATUS='REPLACE', IOSTAT=info)
          if (info /= 0) return
          write(99,*) sOSBW%mCnt,sOSBW%mE0,sOSBW%mDE,sOSBW%mBeta,hwLen
	  if (sCX) then
             write(99,*) HWCX
          else
	     write(99,*) HW
          end if
          close(99)
          saveHWLocal = .true.
     end if

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function loadHWLocal()

     integer ::  info

     loadHWLocal = .false.
     if (srMode) then
          open(99,FILE=fHW, STATUS='OLD', FORM='UNFORMATTED',   &
               IOSTAT=info)
          if (info /= 0) return
          read(99) sOSBW%mCnt,sOSBW%mE0,sOSBW%mDE,sOSBW%mBeta,hwLen
          if (allocHW()) then	      
	      if (sCX) then
	          read(99) HWCX	        
              else
	          read(99) HW                         
              end if             
              loadHWLocal = .true.
           end if
           close(99)
      else
          open(99,FILE=fHW, STATUS='OLD', IOSTAT=info)
          if (info /= 0) return
          read(99,*) sOSBW%mCnt,sOSBW%mE0,sOSBW%mDE,sOSBW%mBeta,hwLen
          if (allocHW()) then	      
	      if (sCX) then
	          read(99,*) HWCX 	          
              else
	          read(99,*) HW                           
              end if              
              loadHWLocal = .true.
           end if
           close(99)
     end if

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
