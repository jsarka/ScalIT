!ccccccccccccccccccccccccccccccccccccccccccccccccc
!c    Subroutine to test Hij used for OSBW       c
!ccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine progHij()

     double precision :: t0, t1
     double precision,allocatable :: E0Ind(:)
     integer, allocatable :: ind(:)

     print *
     print *, '*****************************************'
     print *, '*      Testing Hij used for OSBW        *'
     print *, '*****************************************'
     print *
     print *, ' Center Energy: ', sOSBW%mE0
     print *, ' Expected Energy Window:[', sOSBW%mE0-sOSBW%mDE,    &
                ',',sOSBW%mE0+sOSBW%mDE
     print *, ' Expected # States within Window:',sOSBW%mCnt
     print *, ' Total Number of States:', myLen
     print *

     call CPU_Time(t0)

     if (initOSBW()) then
        allocate(E0Ind(hwlen),ind(hwlen))
        call getOSBWE0Ind(sOSBW%mE0,sOSBW%mDE,hwlen,ind,E0Ind)
        print *, ' Window Size:', hwLen
        print *, ' Indices of Window States:'
        print *, ind
        print *, ' Approximate Energies of Window States:'
        print *, E0Ind
        print *, ' Window Matrix:'
        print *, HW
        deallocate (E0Ind, ind)
        call CPU_Time(t1)
        print *, '  Time to calculate Hij for Window:',t1-t0
     else
        print *, ' Error in calculating Hij for Window:'
     end if
     call finalOSBW()    

     print *, '        Finish Hij Testing '
     print *, '***************************************'
     print *
     
end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine progHijFile(svMode,fname)
     logical,intent(IN) :: svMode
     character(len=*), intent(IN) :: fname

     double precision :: t0, t1
     double precision,allocatable :: E0Ind(:)
     integer, allocatable :: ind(:)
     
     logical :: ltemp, savedata

     print *
     print *, '*****************************************'
     print *, '*      Testing Hij used for OSBW        *'
     print *, '*****************************************'
     print *
     print *, ' Center Energy: ', sOSBW%mE0
     print *, ' Expected Energy Window:[', sOSBW%mE0-sOSBW%mDE,    &
                ',',sOSBW%mE0+sOSBW%mDE
     print *, ' Expected # States within Window:',sOSBW%mCnt
     print *, ' Total Number of States:', myLen
     print *

     call CPU_Time(t0)

     if (initOSBW()) then
        call CPU_Time(t1)
        print *, '  Time to calculate Hij for Window:',t1-t0
        allocate(E0Ind(hwlen),ind(hwlen))
        call getOSBWE0Ind(sOSBW%mE0,sOSBW%mDE,hwlen,ind,E0Ind)
        print *, ' Window Size:', hwLen
        print *, ' Indices of Window States:'
        print *, ind
        print *, ' Approximate Energies of Window States:'
        print *, E0Ind
        print *, ' Window Matrix elements are stored in:', fname
        ltemp = saveData(hwlen**2,HW,svMode,fname)
        IF (.NOT. ltemp) THEN  ! Error Checking -Corey Petty
          PRINT *, "Problem with saving ", fname
        END IF
        deallocate (E0Ind, ind)
     else
        print *, ' Error in calculating Hij for Window:'
     end if
     call finalOSBW()

     call CPU_Time(t1)     

     print *, '   Total Time to calculate/save Hij:', t1-t0
     print *, '        Finish Hij Testing '
     print *, '***************************************'
     print *
     
end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
