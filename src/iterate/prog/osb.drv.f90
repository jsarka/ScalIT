! Driver program for Eigenvalue problem of Tri-atomic molecule
!
program prog_diag
    use osbw 
    
    print *, '**********************************************'
    print *, '     EigValues  of  Molecules  Using PIST     '
    print *, '**********************************************'
    print *
    print *, ' Get Input parameters from STDIO'
    
    call readOSB()

    print *, 'Finish Read Input parameter'

    if (initOSB()) THEN
        call printOSB()
         print *, '  Loading Initial Data .......'
        if (loadInitData()) then
            print *, '  Call ProgDiag subroutine ........'   
            if (progDiag()) then

               print *, '  Get eigenvalues via PIST ..........'
               call progPistEiG()
            else
               print *, '  Error in Block Jacobi Diagonalization!'
            end if  
        else
            print *, '  Error in loading initial data!'
        end if

    else
        print *, '  Memory Allocation error!'
    end if
    
    call finalOSB()

end
!.......................................................................
