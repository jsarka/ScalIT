!
! Driver program for Eigenvalue problem of Tri-atomic molecule
!

program prog_diag
    use mosbp 
    implicit none    
    include 'mtan.dat.h'
    integer :: i, ierr
    logical :: log0

    if (initMOSB()) THEN

       if (id==rootID) then
           print *, '**********************************************'
           print *, '     EigValues  of  Molecules  Using PIST     '
           print *, '**********************************************'
           print *
           print *, ' Get Input parameters from STDIO'

       !    call printMOSB()
           print *, '  Loading Initial Data .......'
        end if

        if (.true.) then
!        if (loadInitDataSeq()) then
        if (loadInitData()) then
            if (id==rootID)   &
                 print *, '  Call ProgDiag subroutine ........'   
            if (progDiag()) then
                if (id==rootID) print *, "eig is stored in ", EIG_FILE
                log0 = saveSeq(.false.,EIG_FILE,eig0)
            else
                print *, '  Error in Block Jacobi Diagonalization!'
            end if  
        else
            print *, '  Error in loading initial data at id=', id
        end if
        end if
    else
        print *, '  MPI Initialization and Memory Allocation error!'
    end if
    
    call finalMOSB()

end

