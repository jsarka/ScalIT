!
! Driver program for Eigenvalue problem of Tri-atomic molecule
!

program prog_diag
    use mosbp 
    implicit none    
    integer :: i, ierr

    if (initMOSB()) THEN

       if (id==rootID) then
           print *, '**********************************************'
           print *, '     EigValues  of  Molecules  Using PIST     '
           print *, '**********************************************'
           print *
           print *, ' Get Input parameters from STDIO'

!           call printMOSB()
           print *, '  Loading Initial Data .......'
        end if

        if (.true.) then
           if (loadInitData()) then   ! Parallel IO
              if (id==rootID)   &
                 print *, '  Call ProgDiag subroutine ........'   
              if (progDiag()) then
                if (id==rootID)  &
                   print *, '  call progPIST! Get eigenvalues via PIST ..........'
                call progPistEig()
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

