!cccccccccccccccccccccccccccccccc
!c    Read input parameters     c
!cccccccccccccccccccccccccccccccc
subroutine readNDPR()
    
   call myreadNDPR(5)

end subroutine

!************************************************
subroutine readNDPRFile(filename)
   character(len=*), intent(IN) :: filename

   open(99, File=Filename, status='OLD')

   call myreadNDPR(99)    

   close(99)

end subroutine

!*************************************************
subroutine myreadNDPR( fd)
   integer, intent(IN) :: fd
   
   double precision :: x1, x2
   
   read(fd, *) BJMax, nGL, nGC
   read(fd, *) lmax, nMax, BJTol, Ecut
   read(fd, *) mass, x1,   x2
   rmin = min(x1,x2); rmax = max(x1,x2)

   read(fd, *) saveMode, useSP
   read(fd, '(A)') xFile    ! store final XYZ DVR points
   read(fd, '(A)') hFile    ! store final H for DVR

   if (useSP)   read(fd, '(A)') inFile 
   
end subroutine

!****************************************************
subroutine printParam()
   integer :: i

   print *
   print *, '================================================================='
   print *, '             Parameters for Non-Direct-Product DVR '
   print *, '================================================================='
   print * 
   print *,' Block Jacobi: Max Tol:', BJTol, ' Max. Iter:',BJMax
   print *,' # Integral pts: Gauss_Legendre:',nGL, '. Gauss_Chebyshev:',nGC
   print *,' j(l) Max:', lmax, '.  # of Orginal DVR pts:',nMax
   print *,' Ecut Energy:', Ecut, '.  Mass:', mass
   print *, ' DVR range:[', rmin, rmax,']'

   print *
   print *, ' Total # of DVR pts:', nDVR, ' Total # of R DVR pts:', nRDVR

   print *
   print *, '================================================================'
   print *, '       l        | # R DVR pts  |  # f(R)*Y(lm) |  Base Address'
   do i = 1, lmax+1
     if (nSize(i)==0) exit
     print *,i-1, nSize(i), nTSize(i), nBase(i)
   end do
   print *, '================================================================='
   print *

end subroutine


!****************************************************
subroutine readVR(filename)
   character(len=*), intent(IN) :: fileName 
   
   double precision :: V0(nMax)
   integer :: i, info

   open(99, file=fileName, status='old')
   read(99,*) spNum

   allocate(spR(spNum),spV(spNum),spM(spNum),stat=info)
   if (info /= 0) then
       close(99); return
   end if

   do i = 1, spNum
         read(99,*) spR(i), spV(i)
   end do
   close(99)

   spyp1 = 0.2D33; spyp2=0.2D33 

   call spline(spNum,spR,spV,spyp1,spyp2,spM)

   call splint(spNum, spR, spV, spM, nMax, X0, V0)

   call myInitVR(V0)

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!**************************************************************
subroutine initVR(fitV)
   external :: fitV

   double precision :: V0(nMax)

   call fitV(nMax, X0, V0)

   call myInitVR(V0)

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!*************************************************************
subroutine myInitVR(V0)
   double precision :: V0(nMax)

   integer :: i

   do i = 1, nMax
        H0(i,i)=H0(i,i) + V0(i)
   end do

   call getRTSize(nSize, nTSize)
   nRDVR = SUM(nSize(1:(lmax+1)))      ! total # of R DVR points
   nDvr  = SUM(nTSize(1:(lmax+1)))     ! total number of basis functions
   nBase(1)=1

   do i=1, lmax
     nBase(i+1)=nBase(i)+nSize(i)
   end do

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
