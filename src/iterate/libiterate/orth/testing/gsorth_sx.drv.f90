!
! Subroutine to test GS-Orthogonalization, Complex version
!

program testGSS
     IMPLICIT NONE
     integer, parameter :: N1 = 1000   ! # of row
     integer, parameter :: N2 = 100   ! # of column
     DOUBLE PRECISION, PARAMETER :: EPSI = 1.0D-14
 
     double complex, dimension (N1, N2) :: A
     double complex, dimension (N2, N2) :: B
     double complex   :: DOT_CX
     double precision :: maxErr, MAXMATI_SX 
     double precision :: tmp1   
     logical :: isOrth_SX
     integer :: i, j

!*************************
     call randMat_CX(N1, N2, A)
   
     print *
     print *, '*************************************'
     print *, ' Testing Complex Version of GS'
     print *, 'Initial Matrix:', N1,'x', N2
!     print *, A(:,1)
!     print *
!     print *, A(:,2)

     CALL GS_FULL_ORTH_SX(N1, N2, A)
     print *
!    print *, 'Matrix after GS ReOrthogonalization:'
!    print *, A(:,1)
!     print *
!     print *, A(:,2)
!     print *

     if (isOrth_SX(N1, N2, A, EPSI)) THEN
        print *, 'A is orthogonal after GS'
     else
        print *, 'Lost orthogonality after GS'
     end if

     call ATA_CX(N1, N2, A, B)
!     print *, 'B=A^H*A:'
!     print *, B
     maxErr = maxMatI_SX(N2, B)        
     print *, 'Maximum error of |A^H*A-I| after GS :', maxErr

!ccccccccccccccccccccccccccccccccccccccccccc
     call randMat_CX(N1, N2, A)

     print *
     print *, ' Testing Complex version of MGS'
!     print *, 'Initial Matrix:', N1, N2
!     print *, A

     CALL MGS_FULL_ORTH_SX(N1, N2, A)
!     print *, 'Matrix after MGS ReOrthogonalization:'
!     print *, A
     
     if (isOrth_SX(N1, N2, A, EPSI)) THEN
        print *, 'A is orthogonal after MGS'
     else
        print *, 'Lost orthogonality after MGS'
     end if

     call ATA_CX(N1, N2, A, B)
     maxErr = maxMatI_SX(N2, B)
     print *, 'Maximum Error of !A^H*A-I| after MGS ', maxErr
  !   print *, 'A^H*A='
  !   print *, B
    print *, '*************************************'
    print *

end




