!cccccccccccccccccccccccccccccccccccccccccc
!c      Do convergence testing            c
!cccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccc
module ndprconv
   IMPLICIT NONE
   character (LEN=128) :: VR_FILENAME
   integer :: lmax, spNum
   integer :: startN, stepN, maxN, num
   double precision :: mass, rmin, rmax
   double precision :: Ecut, maxErr

   double precision, allocatable :: spR(:), spV(:), spM(:)

contains

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!logical function  init()
!   integer :: info
!
!   init = .FALSE.
!   allocate(spR(spNum), spV(spNum), spM(spNum), stat=info)   
!   if (info /= 0) return
!   init = .TRUE.
!
!end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine final()
   if (allocated(spR))  deallocate(spR)
   if (allocated(spV))  deallocate(spV)
   if (allocated(spM))  deallocate(spM)
end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

logical function readVR(filename)
   character(len=*), intent(IN) :: fileName 
   integer :: i, info 
   double precision, parameter :: spyp1 = 3.0D35, spyp2=3.0D35

   readVR = .false.
   open(99, file=filename, status='old')
   read(99,*) spNum

   allocate(spR(spNum),spV(spNum),spM(spNum),stat=info)
   if (info /= 0) then
       close(99); return
   end if

   do i = 1, spNum
         read(99,*) spR(i), spV(i)
   end do
   close(99)

   call spline(spNum,spR,spV,spyp1,spyp2,spM)
   readVR = .true.

end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine calVR(N, X0, V0)
   integer, intent(IN) :: N
   double precision, intent(IN) :: X0(N)
   double precision, intent(OUT):: V0(N)

   call splint(spNum, spR, spV, spM, N, X0, V0)
end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine readInput()
    
   call myread(5)

end subroutine

!**********************************************************
subroutine readFile(filename)
   character(len=*), intent(IN) :: filename

   open(99, File=Filename, status='OLD')
   call myread(99)    
   close(99)

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine myRead(fd)
   integer, intent(IN) :: fd

   read(fd, *) lmax, num, maxErr
   read(fd, *) rmin, rmax, mass, Ecut
   read(fd, *) startN, stepN, maxN
   read(fd, '(A)') VR_FILENAME
  
end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine printParam()

   print *, ' lmax:', lmax, 'Ecut:', Ecut
   print *, ' num:', num, ' Max. Conv Error:', maxErr
   print *, ' Mass:', mass, ' Rmin:', rmin, ' Rmax:', rmax
   print *, ' StartN:', startN, ' StepN:', stepN, ' MaxN:', maxN
   print *, ' # of initial 1D potential: ', spNum

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function doConv()
    double precision :: ER0(num), ER1(num), err, getMax
    integer :: N, info, i
    double precision, allocatable :: X0(:), V0(:), H0(:,:)

    doConv = 0
    do N=startN, maxN, stepN
         allocate(X0(N), V0(N), H0(N,N), stat=info)
         if (info/=0) return
         call Sinc2_XH(N, rMin, rMax, mass, X0, H0)
         call calVR(N, X0, V0)
         do i = 1, N
             H0(i,i)=H0(i,i)+V0(i)
         end do
         if (.NOT. getER0(N, X0, H0, ER1))  return   
         print *, 'N=', N, ' Eig:'
         print *, ER1(1:num)

         if (N==startN) then
            ER0(1:num)=ER1(1:num)
            err = maxErr + 1.0D0
         else
            err = getMax(num, ABS(ER1-ER0))
            print *, ' N=', N, ' Error:', err
            if (err < maxErr) exit
         end if
         deallocate(X0, V0, H0)
    end do 

    if (allocated(X0)) deallocate(X0)
    if (allocated(H0)) deallocate(H0)
    if (allocated(V0)) deallocate(V0)

    if (err < maxErr) then
        doConv = N
    else
        doConv = -N+stepN
    end if

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

logical function getER0(N, X0, H0, ER0)
    integer, intent(IN) :: N
    double precision, intent(IN)  :: X0(N), H0(N, N)
    double precision, intent(OUT) :: ER0(num)
    
    integer, parameter :: WSCALE = 6

    double precision :: hMat(N,N),dHMat(N), work(WSCALE*N)
    integer :: j, i, lwork, info, addr, nSize
    double precision, allocatable :: EH0(:)
    integer, allocatable :: ind(:)    

    lwork = WSCALE*N;    getER0 = .false.

    hMat(1:N,1:N) = H0(1:N, 1:N)
    call DSYEV('N', 'U', N, Hmat, N, DHMAT, work, lwork, info)
    if (info /= 0) return
    nSize = 0
    do i=1,N 
       if (dhMat(i) <= Ecut)  nSize = nSize + 1
    end do
    if (nSize == 0) return

    allocate(EH0((lmax+1)*nSize), stat=info)
    if (info /= 0) return

    addr = 0
    do j=0, lmax
      hMat(1:N,1:N) = H0(1:N, 1:N)
      do i=1, N
         hMat(i,i) = Hmat(i,i) + 0.5D0*j*(j+1)/(Mass*X0(i)*X0(i))
       end do

       call DSYEV('N', 'U', N, hMat, N, DHMAT, work, lwork, info)
       if (info /= 0) return
       nSize = 0
       do i=1, N
          if (dhMat(i) <= Ecut) then
             addr = addr + 1
             EH0(addr)=DHMAT(i)
             nSize = nSize + 1
          end if
       end do
       if (nSize==0) exit
    end do
 
    allocate(ind(addr), stat=info)
    if (info/=0) then
       deallocate(EH0); return
    end if

    call AscOrder(addr, EH0(1:addr), ind)
    addr = min(addr, num)
    do i = 1, addr
       ER0(i) = EH0(ind(i))
    end do
    getER0 = .TRUE.

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

end module
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
