!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c            IO routines to calculate wave3d                  c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine myreadW3(fd)
   integer, intent(IN) :: fd
 
   integer :: getPjmSize, get3Size
   integer :: i,  ind
   double precision :: tmp
 
   read(fd, *)  JTol, parity
   read(fd, *)  Mass(1:3)
   read(fd, *) NState, gType, sType, kNum
   if (NState>0) then     
      do ind = 1, NState
         NSInd(ind)=ind
      end do
   else
      NState=-NState
      read(fd,*) (NSInd(ind), ind=1,NState)
   end if
   read(fd,*) fOut

   do i = 1, 3
      read(fd, *) NR(i), RanMin(i), RanMax(i)
      if (RanMin(i)>RanMax(i)) then
         tmp=RanMin(i); RanMin(i)=RanMax(i); RanMax(i)=tmp
      end if
   end do

   ind = 3
   if ( .NOT. gType) then          ! using (theta)
      if ((RanMin(ind)<Theta0(1)) .OR. (RanMin(ind)>Theta0(2)))  &
            RanMin(ind) = Theta0(1)
      if ((RanMax(ind)>Theta0(2)) .OR. (RanMax(ind)<Theta0(1)))  &
            RanMax(ind) = Theta0(2) 
   else                        ! using  cos(theta)
      if ((RanMin(ind)<XYR0(1)) .OR. (RanMin(ind)>XYR0(2)))  &
            RanMin(ind) = XYR0(1)
      if ((RanMax(ind)>XYR0(2)) .OR. (RanMax(ind)<XYR0(1)))  &
            RanMax(ind) = XYR0(2) 
   end if

   read(fd,*) NMax(1), NS(1), fVlr 
   read(fd,*) NMax(2), NS(2), fVBr 
   read(fd,*) jmax,    NS(3), fVTh 
   read(fd,*) fVP

   p1Nmax = getPjmSize(Jmax)
   NMax(3) = get3Size(parity, JTol, JMax)
   if ((NS(3)<1).OR.NS(3)>NMax(3)) NS(3)=NMax(3)

   if (JTol>jMax) then
     mMax = jMax
   else
     mMax = JTol
   end if

   nTotal=NS(1)*NS(2)*NS(3)

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine printW3()
             
   integer :: i

   write (*,*)
   write(*, *) ' cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
   write(*, *) ' c   Parameter for Triatom Molecule in Jacobi Coordinator   c'
   write(*, *) ' cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc' 
   write(*, 10)  jtol, parity
   write(*, 20)  Mass(1), Mass(2), Mass(3)
   write(*, *)  ' Number of Interested States:', NState
   write(*, *)  ' Indices for Interested States:'
   write(*, *)  (NSInd(I), I=1, NState)
   write(*, 30) '  File to save wave function: ', fOut

   print *, ' Grid information for the parameters:'
   write(*, 40) '  lr:   ', NR(1),RanMin(1), RanMax(1)
   write(*, 40) '  Br:   ', NR(2),RanMin(2), RanMax(2)
   write(*, 40) '  Angle:', NR(3),RanMin(3), RanMax(3)

   if (gType) then
       write(*, *) ' Using cos(theta) for degree of Angle!'
       write(*, 40) '      lr:   ', NR(1),RanMin(1), RanMax(1)
       write(*, 40) '      Br:   ', NR(2),RanMin(2), RanMax(2)
       write(*, 40) ' cos(theta):', NR(3),RanMin(3), RanMax(3)
   else
       write(*, *) ' Using Degree theta for degree of Angle!' 
       write(*, 40) '   lr: ', NR(1),RanMin(1), RanMax(1)
       write(*, 40) '   Br: ', NR(2),RanMin(2), RanMax(2)
       write(*, 40) ' Theta:', NR(3),RanMin(3), RanMax(3)
   end if

   print *
   print *, ' File Names to store original Eigenvectors'

   write(*, 50) '  Eigenvectors for lr: ', NS(1), fVlr
   write(*, 50) '  Eigenvectors for Br: ', NS(2), fVBr
   if (NS(3)/=NMax(3)) then
       write(*,50) ' Eienvectors for Angle:',NS(3),fVTH
   else
	print *, ' No contraction for the Eienvectors of Angle:'
   end if
   write(*, 50) '  Eigenvectors for PIST:', Ntotal, fVP
   write(*, 60)  jmax, p1NMax, NMax(3)

   10 FORMAT('  JTol:',I5,2x,'Parity:', L7)
   20 FORMAT('  Mass:A=',F15.7, '. B=',F15.7, '. C=',F15.7)
   30 FORMAT(A,A)
   40 FORMAT(A, ' Grid Size:',I6, '.  Grid Range:[',F9.4,',',F9.4,']')
   50 FORMAT(A, ' Size:',I10, '. File Name: ', A)
   60 FORMAT('  Jmax:', I10, '.  # of (jk):',I10,'.  Size of Angle Basis:',I6)

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
