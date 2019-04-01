!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Subroutine to read/get eigenvectors for wavefunction  c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Get the wavefunction from Sinc2_DVR basis          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function getPartDVRV_Sinc2(filename, Nmax, NS, Rmin, Rmax, VR)
   implicit none
   character(len=*), intent(IN) :: filename
   integer, intent(IN)  :: NMax, NS
   double precision, intent(OUT) :: Rmin, Rmax, VR(Nmax, NS)

   double precision :: mass
   integer ::  info, N0, i, j, lwork
   double precision, allocatable :: E0(:), V0(:,:),X0(:), X1(:,:), work(:)

   getPartDVRV_Sinc2 = .false.

   if ((Nmax<NS).OR.(Nmax<2)) return
  
   open(99, file=filename, status='old', form='unformatted', iostat=info)
   if (info /= 0) return

   read(99) N0, Rmin, Rmax, mass

   if (N0/=Nmax) then
      close(99); return
   end if

   lwork=3*NS
   allocate(X0(Nmax),E0(Nmax), V0(Nmax,NS),X1(NS,NS), work(lwork), stat=info)
   if (info/=0) then
      close(99); return
   end if

   read(99) X0(1:Nmax), E0(1:Nmax), V0(1:Nmax,1:NS)
   close(99)

   do i = 1, NS
      do j = i, NS
         X1(i,j) = sum(V0(1:NMax,i)*X0(1:NMax)*V0(1:Nmax,j))
         X1(j,i) = X1(i,j)
      end do
   end do

   call  DSYEV('V', 'U', NS, X1, NS, X0, WORK, lwork, info ) 
   call  DGEMM('N','N',NMAX,NS,NS,1.0D0,V0,NMAX,X1,NS,0.0D0,VR,NMax)

   deallocate(X0, E0, X1, V0, work)
   getPartDVRV_Sinc2 = .true.

end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Get the wavefunction from Sinc2_DVR basis          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function getDVRV_Sinc2(filename, NR, R, NS, Vout)
   implicit none
   character(len=*), intent(IN) :: filename
   integer, intent(IN)  :: NR, NS
   double precision, intent(IN)  :: R(NR)
   double precision, intent(OUT) :: Vout(NS, NR)

   double precision :: Rmax, Rmin, dr, mass
   integer :: NMax, info, I, I0, j, lwork
   double precision, allocatable :: E0(:), X0(:),V0(:,:),VR(:,:),X1(:,:),work(:)

   getDVRV_Sinc2 = .false.

   open(99, file=filename, status='old', form='unformatted', iostat=info)
   if (info /= 0) return

   read(99) NMax, Rmin, Rmax, mass

   if ((Nmax<NS).OR. (Nmax)<2) then
      close(99); return
   end if

   lwork=3*NS
   allocate(X0(Nmax),E0(Nmax),V0(Nmax,NS),VR(Nmax,NS),X1(NS,NS),work(lwork),stat=info)
   if (info/=0) then
      close(99); return
   end if

   read(99) X0(1:Nmax), E0(1:Nmax), V0(1:Nmax,1:NS)
   close(99)
  
   do i = 1, NS
      do j = i, NS
         X1(i,j) = sum(V0(1:NMax,i)*X0(1:NMax)*V0(1:Nmax,j))
         X1(j,i) = X1(i,j)
      end do
   end do

   call  DSYEV('V', 'U', NS, X1, NS, X0, WORK, lwork, info ) 
   call  DGEMM('N','N',NMAX,NS,NS,1.0D0,V0,NMAX,X1,NS,0.0D0,VR,NMax)

   dr = (Rmax-Rmin)/Nmax
   do i = 1, NR
      I0 = (R(I)-Rmin)/dr + 1
      if (I0 < 1)    I0 = 1
      if (I0 > Nmax) I0 = Nmax
      Vout(1:NS,I) = VR(I0, 1:NS)
   end do
   
   deallocate(X0, E0, V0, VR, X1, work)
   getDVRV_Sinc2 = .true.

end 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Get the wavefunction of the contracted basis      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function getConV(filename, N, M, V)
   implicit none
   character(len=*), intent(IN) :: filename
   integer, intent(IN) :: N, M
   double precision, intent(OUT) :: V(N, M)

   double precision :: E0(N)
   integer :: info 

   getConV = .false.
   if (N<M) return

   open(99,file=filename,status='old',form='unformatted',iostat=info)
   if (info/=0) return

   read(99) info
   if (info==N) then
      read(99) E0(1:N), V(1:N,1:M)
      getConV = .true.
   end if
   close(99)
end 


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c              Get the PIST wavefunction                       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

logical function getPISTV(filename, N, NS, V)
   implicit none
   character(len=*), intent(IN) :: filename
   integer, intent(IN) :: N, NS
   double precision, intent(OUT) :: V(N, NS)

   double precision, allocatable :: VP(:,:), V0(:,:), E0(:), work(:)
   integer :: N0, NP, info, lwork

   getPISTV = .false.

   open(99,file=filename,status='old',form='unformatted',iostat=info)
   if (info/=0) return

   read(99) NP, N0
   if ((N0/=N).OR.(NP<NS)) then
      close(99); return
   end if

   lwork = 3*NP
   allocate(E0(NP), VP(NP, NP), V0(N,NP), work(lwork), stat=info)
   if (info/=0) then
      close(99); return
   end if 

   read(99) E0(1:NP), VP(1:NP,1:NP)
   read(99) V0(1:N,1:NP)
   close(99)

   call DSYEV('V','U', NP, VP, NP, E0, work, lwork, info)

   if (info==0) then   ! V(N,NS)=V0(N,NP)*VP(NP,NS), NS<NP
       call DGEMM('N','N', N, NS, NP, 1.0D0, V0, N, VP, NP, 0.0D0, V, N)
       print *
       print *, ' ======  PIST Eigen Values ======'
       print *, E0(1:NS)
       print *, ' ================================'
       getPISTV = .true.
   end if

   deallocate(E0, V0, VP, work)

end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c              Get the PIST wavefunction                       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

logical function getPISTV_Index(filename, N, NS, NSind, V)
   implicit none
   character(len=*), intent(IN) :: filename
   integer, intent(IN) :: N, NS, NSInd(NS)
   double precision, intent(OUT) :: V(N, NS)

   double precision, allocatable :: VP(:,:),V0(:,:),V1(:,:),work(:),E0(:)
   integer :: NP, info, N0, i, lwork, ind

   getPISTV_Index = .false.

   open(99,file=filename,status='old',form='unformatted',iostat=info)
   if (info/=0) return

   read(99) NP, N0

   if ((N0/=N).OR.(NP<NS)) then
      close(99); return
   end if

   lwork = 3*NP
   allocate(E0(NP),VP(NP,NP),V0(N,NP),V1(NP,NS),work(lwork),stat=info)
   if (info/=0) then
      close(99); return
   end if

   read(99) E0(1:NP), VP(1:NP,1:NP)
   read(99) V0(1:N,1:NP)
   close(99)

   call DSYEV('V','U', NP, VP, NP, E0, work, lwork, info)

   if (info==0) then
       print *
       print *, ' ======  PIST Eigen Values ======'
       print *, (E0(NSInd(i)),i=1,NS)
       print *, ' ================================'

       do i = 1, NS
         ind = NSInd(i)
         if ( ind < 1 )  ind = 1
         if ( ind > NP ) ind = NP
         V1(1:NP,i) = VP(1:NP, Ind)
       end do 
       
          ! V(N,NS)=V0(N,NP)*V1(NP,NS)
       call DGEMM('N','N', N, NS, NP, 1.0D0, V0, N, V1, NP, 0.0D0, V, N)

     getPISTV_Index = .true.

   end if

   deallocate(E0, V0, VP, V1, work)

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
