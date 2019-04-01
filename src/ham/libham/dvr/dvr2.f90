
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Some utility functions to use DVR                    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical FUNCTION DVR2(N, X10, H1, POT, M, X0, E0, H0)
  implicit none
  integer, intent(in) :: N, M
  double precision, intent(in)    :: X10(N)
  double precision, intent(INOUT) :: H1(N, N)
  external :: POT
  double precision, intent(out)   :: X0(M), E0(M), H0(M,M)

  logical :: DVR_Step2, DVR_Step3
  double precision :: HEig(N)

  DVR2 = .false.
  if ( N < M) return

  DVR2 = DVR_STEP2(N, X10, H1, Pot, HEig)

  if (DVR2) then
     E0(1:M) = HEIG(1:M)
     DVR2 = DVR_STEP3(N, M, X10, H1, E0, X0, H0)
  end if

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical FUNCTION DVR2SF(N, X10, H1, POT, M, X0, E0, H0, fname)
  implicit none
  integer, intent(in) :: N, M
  double precision, intent(in)    :: X10(N)
  double precision, intent(INOUT) :: H1(N, N)
  external :: POT
  double precision, intent(out)   :: X0(M), E0(M), H0(M,M)
  character(len=*), intent(IN)    :: fname

  logical ::  DVR_Step2, DVR_Step3
  double precision :: HEig(N)
  integer :: info

  DVR2SF = .false.
  if ( N < M) return

  DVR2SF = DVR_STEP2(N, X10, H1, Pot, HEig)     

  if (DVR2SF) then
     open(99, file=fname, status='Replace', form='unformatted', iostat=info)     
     if (info/=0) then
            DVR2SF = .false.; return
     end if  
     write(99) N     
     write(99) X10(1:N), HEig(1:N), H1(1:N, 1:N)     
     close(99)

     E0(1:M) = HEIG(1:M)
     DVR2SF = DVR_STEP3(N, M, X10, H1, E0, X0, H0)
  
  end if

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical FUNCTION DVR2RF( M, X0, E0, H0, fname)
  implicit none
  integer, intent(in) :: M
  double precision, intent(out)   :: X0(M), E0(M), H0(M,M)
  character(len=*), intent(IN)    :: fname

  logical :: DVR_Step3
  double precision, allocatable :: H1(:,:), X10(:), HEig(:)
  integer :: info, N 

  DVR2RF = .false.

  open(99, file=fname, status='Old', form='unformatted', iostat=info)
     
  if (info/=0) return
   
  read(99) N  
  allocate(H1(N,M),X10(N),HEig(N), stat=info)     
  if ( (info /= 0) .OR. ( N < M ) ) then
     close(99); return
  end if

  read(99) X10(1:N), HEig(1:N), H1(1:N, 1:M)     
  close(99)
     
  E0(1:M) = HEIG(1:M)
  DVR2RF = DVR_STEP3(N, M, X10, H1, E0, X0, H0)

  deallocate(H1, X10, Heig)

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                  Complex      Version                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical FUNCTION DVRCx2(N, X10, H1, POT, M, X0, E0, H0)
  implicit none
  integer, intent(in) :: N, M
  double precision, intent(in)  :: X10(N)
  double complex, intent(INOUT) :: H1(N, N)
  external :: POT
  double precision, intent(out) :: X0(M), E0(M)
  double complex, intent(OUT)   :: H0(M,M)

  logical :: DVRCx_Step2, DVRCx_Step3
  double precision :: HEig(N)

  DVRCx2 = .false.
  if ( N < M) return

  DVRCx2 = DVRCx_STEP2(N, X10, H1, Pot, HEig)

  if (DVRCx2) then
     E0(1:M) = HEIG(1:M)
     DVRCx2 = DVRCx_STEP3(N, M, X10, H1, E0, X0, H0)
  end if

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical FUNCTION DVRCx2SF(N, X10, H1, POT, M, X0, E0, H0, fname)
  implicit none
  integer, intent(in) :: N, M
  double precision, intent(in)  :: X10(N)
  double complex, intent(INOUT) :: H1(N, N)
  external :: POT
  double precision, intent(out) :: X0(M), E0(M)
  double complex, intent(OUT)   :: H0(M,M)
  character(len=*), intent(IN)  :: fname

  logical ::  DVRCx_Step2, DVRCx_Step3
  double precision :: HEig(N)
  integer :: info

  DVRCx2SF = .false.
  if ( N < M) return

  DVRCx2SF = DVRCx_STEP2(N, X10, H1, Pot, HEig)     

  if (DVRCx2SF) then
     open(99, file=fname, status='Replace', form='unformatted', iostat=info)     
     if (info/=0) then
            DVRCx2SF = .false.; return
     end if  
     write(99) N     
     write(99) X10(1:N), HEig(1:N), H1(1:N, 1:N)     
     close(99)

     E0(1:M) = HEIG(1:M)
     DVRCx2SF = DVRCx_STEP3(N, M, X10, H1, E0, X0, H0)
  
  end if

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical FUNCTION DVRCx2RF( M, X0, E0, H0, fname)
  implicit none
  integer, intent(in) :: M
  double precision, intent(out) :: X0(M), E0(M)
  double complex, intent(OUT)   :: H0(M,M)
  character(len=*), intent(IN)  :: fname

  logical :: DVRCx_Step3
  double complex, allocatable   :: H1(:,:)
  double precision, allocatable :: X10(:), HEig(:)
  integer :: info, N 

  DVRCx2RF = .false.

  open(99, file=fname, status='Old', form='unformatted', iostat=info)
     
  if (info/=0) return
   
  read(99) N  
  allocate(H1(N,M),X10(N),HEig(N), stat=info)     
  if ( (info /= 0) .OR. ( N < M ) ) then
     close(99); return
  end if

  read(99) X10(1:N), HEig(1:N), H1(1:N, 1:M)     
  close(99)
     
  E0(1:M) = HEIG(1:M)
  DVRCx2RF = DVRCx_STEP3(N, M, X10, H1, E0, X0, H0)

  deallocate(H1, X10, Heig)

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
