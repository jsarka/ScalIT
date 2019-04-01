!
! Subroutine to test GS-Orthogonalization in MPI
!

program testGS_MPI
     IMPLICIT NONE
     include 'mpif.h'
     integer, parameter :: N1 = 100   ! # of row
     integer, parameter :: N2 = 50    ! # of column
     DOUBLE PRECISION, PARAMETER :: EPSI = 1.0D-14
 
     double precision, dimension (N1, N2) :: A
     double precision, dimension (N2, N2) :: B
     double precision :: maxErr, MAXMATI    
     logical :: isOrth_MPI
     integer :: myid, ierr, node
     double precision :: mt0, mt1, mt2, ct0, ct1, ct2

!*****************************************

     CALL CPU_Time(ct0)

     CALL MPI_INIT(ierr)
     CALL MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
     CALL MPI_COMM_SIZE(MPI_COMM_WORLD, node, ierr)
     mt0 = MPI_WTime()

     IF (myid == 0) THEN
        print *
        print *, '************************************'
        print *, '     Testing Real Version of GS'
        print *, 'Initial Matrix:', N1*NODE, 'x', N2
   !     print *, A
     END IF

     call random_number(A)
     mt1 = MPI_WTime()
     CALL CPU_Time(ct1)
     CALL GS_FULL_ORTH_MPI(N1, N2, A, IERR)
     mt2 = MPI_WTime()
     CALL CPU_Time(ct2)
!     print *, 'Matrix after GS ReOrthogonalization:'
!     print *, A


     if (isOrth_MPI(N1, N2, A, EPSI)) THEN
        IF (myid == 0) THEN
            print *, 'A is orthogonal after GS'
        END IF
     else
        IF (myid == 0) THEN
            print *, 'Lost orthogonality after GS'
        END IF
     end if

     call DOTPROD_MV_MPI(N1, N2, A, N2, A, MPI_COMM_WORLD, B, ierr)
     maxErr = maxmatI(N2, B)        
     IF (myid == 0) THEN
          print *, 'Maximum error of |A^T*A-I| after GS :', maxErr
          print *, '  Time for Gram Schmidt Orth '
          print *, ' CPU Time:', ct2-ct1, '    MPI WTime:', mt2-mt1
     END IF
!****************************************
     call random_number(A)
     IF (myid == 0) THEN
         print *
         print *, '      Testing Real version of GS'
         print *, 'Initial Matrix:', N1*NODE, N2
    !     print *, A
     END IF

     mt1 = MPI_WTime()
     CALL CPU_Time(ct1)
     CALL MGS_FULL_ORTH_MPI(N1, N2, A, IERR)
     mt2 = MPI_WTime()
     CALL CPU_Time(ct2)
!     print *, 'Matrix after MGS ReOrthogonalization:'
!     print *, A
     
     if (isOrth_MPI(N1, N2, A, EPSI)) THEN
        IF (myid == 0) THEN
           print *, 'A is orthogonal after MGS'
         END IF
     else
        IF (myid == 0) THEN
           print *, 'Lost orthogonality after MGS'
        END IF
     end if

     call DOTPROD_MV_MPI(N1, N2, A, N2, A, MPI_COMM_WORLD, B, ierr)     
     maxErr = maxMatI(N2, B)
     IF (myid == 0) THEN
         print *, 'Maximum Error of |A^T*A-I| after MGS ', maxErr
         print *, '  Time for Modified Gram Schmidt Orth '
         print *, ' CPU Time:', ct2-ct1, '    MPI WTime:', mt2-mt1
         print *, '****************************************'
         mt2 = MPI_WTime()
         CALL CPU_Time(ct2)
         print *, '  Total Time for Gram Schmidt Testing '
         print *, ' CPU Time:', ct2-ct0, '    MPI WTime:', mt2-mt0
         print *, '***********************************************'
         print *
     END IF

     CALL MPI_FINALIZE()
end




