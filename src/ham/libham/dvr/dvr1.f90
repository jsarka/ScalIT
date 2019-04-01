
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Some utility functions to use DVR                    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical FUNCTION DVR1(N, X1, H1, POT, M, X0, E0, H0)
  implicit none
  integer, intent(in) :: N, M
  double precision, intent(inout) :: X1(N,N), H1(N, N)
  external :: POT
  double precision, intent(out)   :: X0(M), E0(M), H0(M,M)

  logical :: DVR_Step1, DVR_Step2, DVR_Step3
  double precision,allocatable :: X10(:), HEig(:)
  integer :: info

  DVR1 = .false.;      if ( N < M) return
  allocate(X10(N),HEig(N), stat=info)
  if (info/=0)  return

  DVR1 = DVR_STEP1(N, X1, H1, X10)

  if (DVR1) then

     DVR1 = DVR_STEP2(N, X10, H1, Pot, HEig)

     if (DVR1) then
         E0(1:M) = HEIG(1:M)
         DVR1 = DVR_STEP3(N, M, X10, H1, E0, X0, H0)
     end if

  end if

  deallocate(X10, HEig)

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical FUNCTION DVR1SF(N, X1, H1, POT, M, X0, E0, H0, fname1, fname2)
  implicit none
  integer, intent(in) :: N, M
  double precision, intent(inout) :: X1(N,N), H1(N, N)
  external :: POT
  double precision, intent(out)   :: X0(M), E0(M), H0(M,M)
  character(len=*), intent(IN)    :: fname1, fname2

  logical :: DVR_Step1, DVR_Step2, DVR_Step3
  double precision, allocatable :: X10(:),HEig(:)
  integer :: info

  DVR1SF = .false.;      if ( N < M) return
  allocate(X10(N),HEig(N), stat=info)
  if (info/=0)  return

  DVR1SF = DVR_STEP1(N, X1, H1, X10)

  if (DVR1SF) then

     open(99, file=fname1, status='Replace', form='unformatted', iostat=info)     
     if (info/=0) then
         DVR1SF = .false.;  deallocate(X10, HEig);  return
     end if
     write(99) N
     write(99) X10(1:N), H1(1:N, 1:N), X1(1:N, 1:N) 
     close(99)

     DVR1SF = DVR_STEP2(N, X10, H1, Pot, HEig)     

     if (DVR1SF) then
        open(99, file=fname2, status='Replace', form='unformatted', iostat=info)     
        if (info/=0) then
            DVR1SF = .false.;  deallocate(X10, HEig); return
        end if  
        write(99) N     
        write(99) X10(1:N), HEig(1:N), H1(1:N, 1:N)     
        close(99)

        E0(1:M) = HEIG(1:M)
        DVR1SF = DVR_STEP3(N, M, X10, H1, E0, X0, H0)
     end if  
  end if

  deallocate(X10, HEig)

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical FUNCTION DVR1RF(POT, M, X0, E0, H0, fname1, fname2)
  implicit none
  integer, intent(in) :: M
  external :: POT
  double precision, intent(out)   :: X0(M), E0(M), H0(M,M)
  character(len=*), intent(IN)    :: fname1, fname2

  logical :: DVR_Step2, DVR_Step3
  double precision, allocatable :: H1(:,:), X10(:), HEig(:)
  integer :: info, N 

  DVR1RF = .false.
   
  open(99, file=fname1, status='Old', form='unformatted', iostat=info)
     
  if (info/=0) return
   
  read(99) N

  allocate(H1(N,N),X10(N),HEig(N), stat=info)
     
  if ((info /= 0) .OR. (N  < M) )then
     close(99); return
  end if

  read(99) X10(1:N), H1(1:N, 1:N)
     
  close(99)
     
  DVR1RF = DVR_STEP2(N, X10, H1, Pot, HEig)

  if (DVR1RF) then
        
     open(99, file=fname2, status='Replace', form='unformatted', iostat=info)
     
     if (info == 0 ) then  
        write(99) N     
        write(99) X10(1:N), HEig(1:N), H1(1:N, 1:N)     
        close(99)

        E0(1:M) = HEIG(1:M)
        DVR1RF = DVR_STEP3(N, M, X10, H1, E0, X0, H0)

     end if
  end if

  deallocate(H1, X10, Heig)

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc





!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c              Complex       Version                          c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical FUNCTION DVRCx1(N, X1, H1, POT, M, X0, E0, H0)
  implicit none
  integer, intent(in) :: N, M
  double complex, intent(inout) :: X1(N,N), H1(N, N)
  external :: POT
  double precision, intent(out) :: X0(M), E0(M)
  double complex, intent(OUT)   :: H0(M,M)

  logical :: DVRCx_Step1, DVRCx_Step2, DVRCx_Step3
  double precision,allocatable :: X10(:), HEig(:)
  integer :: info

  DVRCx1 = .false.;      if ( N < M) return
  allocate(X10(N),HEig(N), stat=info)
  if (info/=0)  return

  DVRCx1 = DVRCx_STEP1(N, X1, H1, X10)

  if (DVRCx1) then

     DVRCx1 = DVRCx_STEP2(N, X10, H1, Pot, HEig)

     if (DVRCx1) then
         E0(1:M) = HEIG(1:M)
         DVRCx1 = DVRCx_STEP3(N, M, X10, H1, E0, X0, H0)
     end if

  end if

  deallocate(X10, HEig)

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical FUNCTION DVRCx1SF(N, X1, H1, POT, M, X0, E0, H0, fname1, fname2)
  implicit none
  integer, intent(in) :: N, M
  double complex, intent(inout) :: X1(N,N), H1(N, N)
  external :: POT
  double precision, intent(out) :: X0(M), E0(M)
  double complex, intent(OUT)   :: H0(M,M)
  character(len=*), intent(IN)  :: fname1, fname2

  logical :: DVRCx_Step1, DVRCx_Step2, DVRCx_Step3
  double precision, allocatable :: X10(:),HEig(:)
  integer :: info

  DVRCx1SF = .false.;      if ( N < M) return
  allocate(X10(N),HEig(N), stat=info)
  if (info/=0)  return

  DVRCx1SF = DVRCx_STEP1(N, X1, H1, X10)

  if (DVRCx1SF) then

     open(99, file=fname1, status='Replace', form='unformatted', iostat=info)     
     if (info/=0) then
         DVRCx1SF = .false.;  deallocate(X10, HEig);  return
     end if
     write(99) N
     write(99) X10(1:N), H1(1:N, 1:N), X1(1:N, 1:N) 
     close(99)

     DVRCx1SF = DVRCx_STEP2(N, X10, H1, Pot, HEig)     

     if (DVRCx1SF) then
        open(99, file=fname2, status='Replace', form='unformatted', iostat=info)     
        if (info/=0) then
            DVRCx1SF = .false.;  deallocate(X10, HEig); return
        end if  
        write(99) N     
        write(99) X10(1:N), HEig(1:N), H1(1:N, 1:N)     
        close(99)

        E0(1:M) = HEIG(1:M)
        DVRCx1SF = DVRCx_STEP3(N, M, X10, H1, E0, X0, H0)
     end if  
  end if

  deallocate(X10, HEig)

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical FUNCTION DVRCx1RF(POT, M, X0, E0, H0, fname1, fname2)
  implicit none
  integer, intent(in) :: M
  external :: POT
  double precision, intent(out)   :: X0(M), E0(M), H0(M,M)
  character(len=*), intent(IN)    :: fname1, fname2

  logical :: DVRCx_Step2, DVRCx_Step3
  double precision, allocatable :: X10(:), HEig(:)
  double complex, allocatable :: H1(:,:)
  integer :: info, N 

  DVRCx1RF = .false.
   
  open(99, file=fname1, status='Old', form='unformatted', iostat=info)
     
  if (info/=0) return
   
  read(99) N

  allocate(H1(N,N),X10(N),HEig(N), stat=info)
     
  if ((info /= 0) .OR. (N  < M) )then
     close(99); return
  end if

  read(99) X10(1:N), H1(1:N, 1:N)
     
  close(99)
     
  DVRCx1RF = DVRCx_STEP2(N, X10, H1, Pot, HEig)

  if (DVRCx1RF) then
        
     open(99, file=fname2, status='Replace', form='unformatted', iostat=info)
     
     if (info == 0 ) then  
        write(99) N     
        write(99) X10(1:N), HEig(1:N), H1(1:N, 1:N)     
        close(99)

        E0(1:M) = HEIG(1:M)
        DVRCx1RF = DVRCx_STEP3(N, M, X10, H1, E0, X0, H0)

     end if
  end if

  deallocate(H1, X10, Heig)

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
