!
! Subroutine to test GS-Orthogonalization
!

program testGS
     IMPLICIT NONE
     integer, parameter :: N1 = 1000   ! # of row
     integer, parameter :: N2 = 500    ! # of column
     DOUBLE PRECISION, PARAMETER :: EPSI = 1.0D-14
 
     double precision, dimension (N1, N2) :: A
     double precision, dimension (N2, N2) :: B
     double precision :: maxErr, MAXMATI    
     logical :: isOrth

!*****************************************
     call random_number(A)
     print *
     print *, '************************************'
     print *, '     Testing Real Version of GS'
     print *, 'Initial Matrix:', N1, 'x', N2
!     print *, A

     CALL GS_FULL_ORTH(N1, N2, A)
!     print *, 'Matrix after GS ReOrthogonalization:'
!     print *, A
     
     if (isOrth(N1, N2, A, EPSI)) THEN
        print *, 'A is orthogonal after GS'
     else
        print *, 'Lost orthogonality after GS'
     end if

     call ATA(N1, N2, A, B)
     maxErr = maxmatI(N2, B)        
     print *, 'Maximum error of |A^T*A-I| after GS :', maxErr

!****************************************
     call random_number(A)
     print *
     print *, '      Testing Real version of GS'
     print *, 'Initial Matrix:', N1, N2
!     print *, A

     CALL MGS_FULL_ORTH(N1, N2, A)
!     print *, 'Matrix after MGS ReOrthogonalization:'
!     print *, A
     
     if (isOrth(N1, N2, A, EPSI)) THEN
        print *, 'A is orthogonal after MGS'
     else
        print *, 'Lost orthogonality after MGS'
     end if

     call ATA(N1, N2, A, B)
     maxErr = maxMatI(N2, B)
     print *, 'Maximum Error of |A^T*A-I| after MGS ', maxErr
     print *, '****************************************'
     print *
end




