!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      Driver program to testing time scalability      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

program prog_test
    use mosbp 
    implicit none    
    integer :: i, ierr
    double precision :: wt1,wt2
    logical :: log0

    if (initMOSB()) THEN
       wt1 = MPI_WTime()
       if (id==rootID) then
           print *, '**********************************************'
           print *, '     MPI_WTime() Tesing of OSB/OSBW program   '
           print *, '**********************************************'
           print *
           print *, ' Get Input parameters from STDIO'

!           call printMOSB()
           print *, '  Loading Initial Data .......'
        end if

        if (loadInitDataSeq()) then   ! without Parallel IO
            if (id==rootID)   &
                 print *, '  Call ProgDiag subroutine ........'   
            log0 = progDiag()

            if (id==rootID)  &
                 print *, '  Call ProgHX subroutine ........'
            call progHX()
           
            if (id==rootID)  &
                 print *, '  Call ProgPX subroutine ........'
            call progPX()

            if (id==rootID)  &
                 print *, '  Call ProgQMR subroutine ........'
            call progQMR()  
        else
            print *, '  Error in loading initial data at id=', id
        end if

        wt2=MPI_WTime()
        if (id==rootID) print *, ' The total MPI WTime for the program:',wt2-wt1        
    else
        print *, '  MPI Initialization and Memory Allocation error!'
    end if

    if (id==rootID) print *, ' ===========  Finish  the program  =============='

    call finalMOSB()

end

