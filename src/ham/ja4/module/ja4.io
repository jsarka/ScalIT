!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Read input parameters either from stdIO or a file  c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine readJA4StdIO()
   call myreadJA4(STDFH)
end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine readJA4File(filename)
   character(len=*), intent(IN) :: filename

   open(99, File=Filename, status='OLD')
   call myreadJA4(99)    
   close(99)

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         Real work in reading input parameters        c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine myreadJA4( fd)
   integer, intent(IN) :: fd

   double precision :: x1, x2

   read(fd, *) JTol, parity 
   read(fd, *) JMAX(1:NA), NGI(1:NA)
   read(fd, *) FcFlag, CbFlag, AbsFlag, useSP, Ecutoff

   read(fd, '(A)') fH0
   read(fd, '(A)') fH0GM

   read(fd, *) MASS(1), RE(1), NDVR(1)
   read(fd, '(A)') fVRlr1

   read(fd, *) MASS(2), RE(2), NDVR(2)
   read(fd, '(A)') fVRlr2

   read(fd, *) MASS(3), RE(3), NDVR(3)
   read(fd, '(A)') fVRBR

   read(fd, *) nDVR(4), ReFlag
   read(fd, '(A)') fRe

   if (useSP) then
      read(fd, '(A)') fSPVRlr1
      read(fd, '(A)') fSPVRlr2
      read(fd, '(A)') fSPVRBR
   end if

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Print out related parameters               c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine printJA4()
   write (*,*)
   write(*, *) '******    Parameters for Jacob-Tetra-Atom Module    *****'
   write(*, *) '                   Grid Information                  '
   write(*, *) '========================================================='
   write(*, *) '        Equ. |    Grid #   |    Mass '
   write(*, *) '========================================================='
   write(*, 10) ' lr1:', RE(1), NDVR(1),  Mass(1)
   write(*, 10) ' lr2:', RE(2), NDVR(2),  MASS(2)
   write(*, 10) ' BR: ', RE(3), NDVR(3),  MASS(3)
   write(*, *) '========================================================='
   write(*, *)
   write(*, 20) jtol, parity
   write(*, 30) jmax(1),jmax(2),jmax(3), jkNum
   write(*, 40) NGI(1), NGI(2), NGI(3)
   write(*, 50) p1NMax, p2NMax, mmNum, mMax
    
   print *
   write (*,5) ' File name to store H0:', fH0
   write (*,5) ' File name to store H0GM:', fH0GM

   write (*,5) ' File name storing diag infor of Vlr1:', fVRlr1
   write (*,5) ' File name storing diag infor of Vlr2:', fVRlr2
   write (*,5) ' File name storing diag infor of VBr:', fVRBr
   write (*,5) ' Save V and E0 of H0RE to binary file:', fRe

   select case (FcFlag)
   case(:FCNONE)
      print *, ' None of radials (r1, r2, R) is fixed.'
   case (FCALL:)
      print *, ' All the radials (r1, r2, R) are fixed.'
   case (FCR1R2)
      print *, ' Radials (r1, r2) are fixed.'
   case (FCBRR1)
      print *, ' Radials (R, r1) are fixed.'
   case (FCBRR2)
      print *, ' Radials (R, r2) are fixed.'
   case (FCBR)
      print *, ' Radial (R) is fixed.'
   case (FCR1)
      print *, ' Radial (R1) is fixed.'
   case (FCR2)
      print *, ' Radial (R2) is fixed.'
   end select

   print *
   select case (CBFlag)
   case(:CBNONE)
      print *, ' None of radial coordinators is combined.'
   case (CBALL:)
      print *, ' All the radial coordinators are combined.'
      print *, ' Cutoff Energy:', Ecutoff
   case (CBR1R2)
      print *, ' Radials (r1, r2) are combined.'
      print *, ' Cutoff Energy:', Ecutoff
   case (CBBRR1)
      print *, ' Radials (R, r1) are combined.'
      print *, ' Cutoff Energy:', Ecutoff
   case (CBBRR2)
      print *, ' Radials (R, r2) are combined.'
      print *, ' Cutoff Energy:', Ecutoff
   end select

   if (useSP) then
       print *, ' Use Spline Function to calculate 1D potential !'
       write (*,5) ' File name for spline function of Vlr1:', fSPVRlr1
       write (*,5) ' File name for spline function of Vlr2:', fSPVRlr2
       write (*,5) ' File name for spline function of VBR:', fSPVRBR
   else
       print *, ' Use Fitting Function to calculate 1D potential !'
   end if

   write(*, *) '=========================================================='
   write(*, *) '********    End Parameters for Jacob-Tetra-Atom    *********'

    5 FORMAT(A, A)
   10 FORMAT(' ', A, F10.5, 2x, I8, 2x, F15.6)
   20 FORMAT(' Total J:',I5,2x,' System Parity:', L3)
   30 FORMAT(' j1max:',I5,2x,' j2max:', I5,2x, 'jmax:',I5, 2x, '# of Angles:',I10)
   40 FORMAT(' Integral parameters: G-Leg1:',I6,2x,' G-Leg2:',I6,2x,' G-Chev Phi:',I6)
   50 FORMAT(' # of Pjm1:',I6,3x,'# of Pjm2:',I6,3x,' # of CG:',I8,3x,' Max|m1-m2|:',I6)

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc





