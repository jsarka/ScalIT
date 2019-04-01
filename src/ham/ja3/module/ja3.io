!cccccccccccccccccccccccccccccccccccccccc
!c          Read input parameters       c
!cccccccccccccccccccccccccccccccccccccccc
subroutine readJA3StdIO()

   call myreadJA3(STDIN)

end subroutine

!************************************************
subroutine readJA3File(filename)
   character(len=*), intent(IN) :: filename

   open(99, File=Filename, status='OLD')
   call myreadJA3(99)    
   close(99)

end subroutine

!*************************************************
subroutine myreadJA3(fd)
   integer, intent(IN) :: fd
   integer :: tmp

   read(fd, *) JTol, parity, JMAX(1), NGI(1)
   read(fd, *) FcFlag, CbFlag, AbsFlag, useSP, Ecutoff

   read(fd, '(A)') fH0
   read(fd, '(A)') fH0GM

   read(fd, *) MASS(1), RE(1), NDVR(1)
   read(fd, '(A)') fVRlr

   read(fd, *) MASS(2), RE(2), NDVR(2)
   read(fd, '(A)') fVRBR

   read(fd, *) nDVR(3), ReFlag
   read(fd, '(A)') fRE

   if (useSP) then
      read(fd, '(A)') fSpVRlr
      read(fd, '(A)') fSpVrBR
   end if

   en(1:2)=0;        A0(1:2)=0.0D0; 
   Rabs0(1:2)=0.0D0; Rabs1(1:2)=0.0D0
   select case (AbsFlag)
   case (ABS_ONE)   ! one absorption potential
        read(fd,*) en(1),A0(1),Rabs0(1),Rabs1(1)
        read(fd,'(A)') fABS
        if (Rabs0(1)>Rabs1(1)) then
            tmp=Rabs0(1); Rabs0(1)=Rabs1(1); Rabs1(1)=tmp
        endif
   case (ABS_TWO)   ! two absorption potentials
        read(fd,*) en(1),A0(1),Rabs0(1),Rabs1(1)
        read(fd,*) en(2),A0(2),Rabs0(2),Rabs1(2)
        read(fd,'(A)') fABS
        if (Rabs0(1)>Rabs1(1)) then
            tmp=Rabs0(1); Rabs0(1)=Rabs1(1); Rabs1(1)=tmp
        endif
        if (Rabs0(2)>Rabs1(2)) then
            tmp=Rabs0(2); Rabs0(2)=Rabs1(2); Rabs1(2)=tmp
        endif
   end select

end subroutine

!****************************************************

subroutine printJA3()
             
   write (*,*)
   write(*, *) '******  Parameters for Jacob-Triatom Module  *****'
  
   print *
   write(*, *) '               Grid Information              '
   write(*, *) '============================================='
   write(*, *) '     Equilli.  |    Grid #   |     Mass   '
   write(*, *) '============================================='
   write(*, 10) ' lr:', RE(1),  NDVR(1), Mass(1)
   write(*, 10) ' BR:', RE(2),  NDVR(2), MASS(2)
   write(*, *) '==============================================='
   write(*, 20) jtol,   parity, jmax(1) 
   write(*, 30) NGI(1), p1NMax
   write(*, 40) jkNum
   print *
   write (*,5) ' Filename to store H0:', fH0
   write (*,5) ' Filename to store H0GM:', fH0GM

   write (*,5) ' Filename for diag info of Vlr:', fVRlr
   write (*,5) ' Filename for diag info of VBr:', fVRBr

   if (ReFlag>0) then
      write (*,5) ' Save V and E0 of H0RE to binary file:', fRE
   else
      if (ReFlag<0) write (*,5) ' Load V and E0 of H0RE from binary file:', fRE
   end if

   select case (FcFlag)
   case(FCNONE)
      print *, ' None of radials (lr, BR) is fixed.'
   case (FCALL:)
      print *, ' All the radials (lr, BR) are fixed.'
   case (FCLR)
      print *, ' Radials (lr) is fixed.'
   case (FCBR)
      print *, ' Radials (R) is fixed.'
   case default
      print *, ' None of radials (lr, BR) is fixed.'
   end select

   print *
   select case (CBFlag)
   case(CBNONE)
      print *, ' None of radial coordinators is combined.'
   case (CBALL)
      print *, ' All the radial coordinators (lr, R) are combined.'
      print *, ' Cutoff Energy:', Ecutoff
   case default
      print *, ' None of radial coordinators is combined.'
   end select

   print *
   if (useSP) then
      print *, ' Using Spline functions to calculate 1D potential'
      write(*,5) ' Filename of Spline VR potential for lr:', fSpVrlr
      write(*,5) ' Filename of Spline VR potentila for BR:', fSpVrBR
   else
      print *, ' Using Fitting functions to calculate 1D potential'
   end if

   print *
   select case (AbsFlag)
   case(ABS_ONE)
       print *, ' Calculate absorption potential(ABS:-W): H = H0 + i(-W) '
       write(*,5) ' Filename to store ABS (-W):', fABS
       print *, ' Polynomials function of ABS for the molecule(W(R)):'
       print *, ' W = A0 * [(R-R0)/(R1-R0)]^n when R > R0; R1>R0'
       print *, ' W = 0                       when R < R0'
       write(*,50) en(1), A0(1), Rabs0(1), Rabs1(1)
   case(ABS_TWO)
       print *, ' Calculate absorption potential(ABS:-W): H = H0 + i(-W) '
       write(*,5) ' Filename to store ABS (-W):', fABS
       print *, ' Polynomials function of reactant ABS(W(R)):'
       print *, ' W = A0 * [(R0-R)/(R1-R0)]^n when R < R0; R1<R0'
       print *, ' W = 0                       when R > R0'
       write(*,50) en(1), A0(1), Rabs0(1), Rabs1(1)
       print *, ' Polynomials function of product ABS(W(R)):'
       print *, ' W = A0 * [(R-R0)/(R1-R0)]^n when R > R0; R1>R0'
       print *, ' W = 0                       when R < R0'
       write(*,50) en(2), A0(2), Rabs0(2), Rabs1(2)
   case default
       print *, ' Do not calculate the absorption potential.'
   end select

   print *
   write(*, *) '==============================================='
   write(*, *) '****    End Parameters for Jacob-Triatom    ***'
   print *
    
   5  FORMAT(' ', A, A)    
   10 FORMAT(' ', A, F10.6, 2x, I8, 4x, F15.9)
   20 FORMAT('  JTol:',I5,4x,'Parity:', L3,4x, 'jmax:', I10)
   30 FORMAT('  # of G-L Integral:',I5,4x, '. # Asso. Legendre:',I10,'.')
   40 FORMAT('  # of (jk) indices:', I10)
   50 FORMAT(' n=', I5, ' A0=',F15.9,' R0=',F15.9,' R1=',F15.9)
end subroutine
!********************************************************
