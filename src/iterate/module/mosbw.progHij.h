!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c          Subroutine to test Hij used for OSBW         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine progHij()

     double precision :: t0, t1
     double precision,allocatable :: E0Ind(:)

     if (id==rootID) then
        print *
        print *, '*****************************************'
        print *, '*      Testing Hij used for OSBW        *'
        print *, '*****************************************'
        t0 = MPI_WTime()
     end if

     if (initMOSBW()) then
        allocate(E0Ind(hwlen))
        call getOSBWE0Ind(E0Ind)

        if (id==rootID) then

   	   call printOSBWInd()

           print *, ' Approximate Energies of Window States:'
           print *, E0Ind

           print *, ' Window Matrix:'
           print *, HW

           t1 = MPI_WTime()
           print *, '  Time to calculate Hij for Window:',t1-t0

        end if

        deallocate (E0Ind)

     else
        if (id==rootID) print *, ' Error in calculating Hij for Window:'
     end if

     call finalMOSBW()    

     if (id==rootID) then
        print *, '        Finish Hij Testing '
        print *, '***************************************'
        print *
     end if
     
end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine progHijFile(svMode,fname)
     logical,intent(IN) :: svMode
     character(len=*), intent(IN) :: fname

     double precision :: t0, t1
     double precision,allocatable :: E0Ind(:)

     if (id==rootID) then
        print *
        print *, '*****************************************'
        print *, '*      Testing Hij used for OSBW        *'
        print *, '*****************************************'

        t0 = MPI_WTime()
     end if


     if (initMOSBW()) then
        allocate(E0Ind(hwlen))
        call getOSBWE0Ind(E0Ind)

        if (id==rootID) then
   	   call printOSBWInd()
          
           print *, ' Approximate Energies of Window States:'
           print *, E0Ind

           print *, ' Window Matrix elements are stored in:', fname
           call saveData(hwlen*hwlen,HW,svMode,fname)

           t1 = MPI_WTime()
           print *, '  Time to calculate Hij for Window:',t1-t0
        end if
        deallocate (E0Ind)
     else
        if (id==rootID) print *, ' Error in calculating Hij for Window:'
     end if

     call finalMOSBW()    

     if (id==rootID) then
        print *, '        Finish Hij Testing '
        print *, '***************************************'
        print *
     end if
     
end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine printOSBWInd()
    integer :: i, j

    print *
    print *, 'ccccccccccccccccccccccccccccccccccccccccccccc'
    print *, 'c     Indices for OSBW Window States        c'
    print *, 'ccccccccccccccccccccccccccccccccccccccccccccc'
    print *
    print *, ' Center Energy: ', sOSBW%mE0
    print *, ' Expected Energy Window:[', sOSBW%mE0-sOSBW%mDE,    &
             ',',sOSBW%mE0+sOSBW%mDE
    print *, ' Expected # States within Window:',sOSBW%mCnt
    print *, ' Total Number of States:', myconf%gN
    print *
    print *, ' # of Nodes:', nNodes, '  ID of Current Node:', id
    print *, ' Window Size:', hwLen, ' Local Window Size:',phwlen
    print *, ' Global of Local Indices of Window States:'
    print *, gInd
    print *, ' Global of Global Indices of Window States:'
    print *, gKIndex
    print *, ' Local of Local Indices of Window States:'
    print *, pInd 
    print *, ' Local of Global Indices of Window States:'
    print *, gPIndex
    print *, ' Node Indices of Window States:'
    print *, nodInd
    print *, ' Node-Count Indices of Window States:'
    print *, gCnt
    print *

    print *, ' Indices for each Windows States'
    do i = 1, hwLen
       print *
       print *, '  Ind    Glb_Ind    level   NodeInd  Grp_Ind   root_Ind  node_Num blkInd sNInd'
       do j = 1, sF
          write(*,10) i, gKIndex(i), j, nodInd(i), grpInd(j,i),rootInd(j,i), &
                      nodeNum(j,i), blkInd(j,i), sNInd(j,i)
       end do
    end do

    10 format (1x,I4,1X,I10,7(1X,I8))
end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
