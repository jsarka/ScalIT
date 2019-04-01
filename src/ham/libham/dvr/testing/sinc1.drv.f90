!c
!c  Testing program for sinc1_dvr
!c  Functions: getsinc_podvr, podvr have been tested.
!c

program sinc1_podvr_test
  implicit none
  double precision, parameter :: MASS_FACTOR=1822.889D0
  integer, parameter :: N = 253
  integer, parameter :: M = 40
  double precision, parameter :: x = 38.0D0
  double precision, parameter :: mass = 10.09D0*MASS_FACTOR
  character(len=16) :: fname='t1.dat' 
  logical, parameter :: isCMU=.true.
 
  double precision :: X0(M), H0(M,M), E0(M), Y0(M), G0(M,M),F0(M)
  
  external :: pot

  logical :: getDVR_sinc1, getDVRSF_Sinc1, getDVRRF_Sinc1

  print *, 'Grid Points:', 2*N+1
  print *, 'Final Grid Points:', M
  print *, 'Coordinate Range:[',-x, x,']'
  print *, 'Mass:', mass

  print *
  print *, ' Testing getDVR_Sinc1 .......'  
  if( getDVR_Sinc1(isCMU,N, X, mass, POT, M, X0, E0, H0) ) THEN
     print *, ' DVR X0:'
     print *,  X0
  !   print *, ' Eig0:', E0
  !   print *, ' Hmatrix:', H0
  else
     print *, ' Error in getDVR_Sinc1'
  end if

  print *
  print *, ' Testing getDVRSF_Sinc1 .......'
  if( getDVRSF_Sinc1(isCMU,N, X, mass, POT, M, Y0, F0, G0, fname) ) THEN
     print *, ' Diff DVR X0:', X0-Y0
 !    print *, ' Diff Eig0:', E0-F0
 !    print *, ' Diff Hmatrix:', H0-G0
  else
     print *, ' Error in getDVRSF_Sinc1'
  end if

  Y0=0;F0=0;G0=0
  print *
  print *, ' Testing getDVRRF_Sinc1 .......'
  if( getDVRRF_Sinc1( M, Y0, F0, G0,fname) ) THEN
     print *, ' Diff DVR X0:', X0-Y0
 !    print *, ' Diff Eig0:', E0-F0
 !    print *, ' Diff Hmatrix:', H0-G0
  else
     print *, ' Error in getDVRRF_Sinc1'
  end if


end 


subroutine POT(N, X, RES)
  implicit none
  double precision, parameter :: E0 = 9.0*24.743267D0  ! in cm^-1
  double precision, parameter :: A0 = 1.0D5              ! in a0

  integer, intent(IN)  :: N
  double precision, intent(in)  :: X(N)
  double precision, intent(out) :: RES(N) 
 
  integer :: i
  double precision :: x6

  do i = 1, N
     if (x(i) == 0.0D0) then
        RES(i) =  - E0
      else
         x6 = X(i)**6     
         RES(i) =  E0*(exp(-A0/X6) - 1)
       end if
  END DO

end
