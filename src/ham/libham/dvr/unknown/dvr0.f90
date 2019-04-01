
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Some utility functions to use DVR                    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical FUNCTION DVR0(N, X1, H1, POT, M, X0, E0, H0)
  implicit none
  integer, intent(in) :: N, M
  double precision, intent(inout) :: X1(N,N), H1(N)
  external :: POT
  double precision, intent(out)   :: X0(M), E0(M), H0(M,M)

  logical :: DVR_Step0, DVR_Step2, DVR_Step3
  double precision,allocatable :: X10(:), HEig(:), H2(:,:)
  integer :: info

  DVR0 = .false.;      if ( N < M) return
  allocate(X10(N),HEig(N),H2(N,N), stat=info)
  if (info/=0)  return

  DVR0 = DVR_STEP0(N, X1, H1, X10, H2)

  if (DVR0) then

     DVR0 = DVR_STEP2(N, X10, H2, Pot, HEig)

     if (DVR0) then
         E0(1:M) = HEIG(1:M)
         DVR0 = DVR_STEP3(N, M, X10, H2, E0, X0, H0)
     end if
  end if

  deallocate(X10, HEig, H2) 

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical FUNCTION DVR0SF(N, X1, H1, POT, M, X0, E0, H0, fname1, fname2)
  implicit none
  integer, intent(in) :: N, M
  double precision, intent(inout) :: X1(N,N), H1(N)
  external :: POT
  double precision, intent(out)   :: X0(M), E0(M), H0(M,M)
  character(len=*), intent(IN)    :: fname1, fname2

  logical :: DVR_Step0, DVR_Step2, DVR_Step3
  double precision,allocatable :: X10(:), HEig(:), H2(:,:) 
  integer :: info

  DVR0SF = .false.;    if ( N < M) return
  allocate(X10(N),HEig(N),H2(N,N), stat=info)
  if (info/=0)  return

  DVR0SF = DVR_STEP0(N, X1, H1, X10, H2)

  if (DVR0SF) then
     open(99, file=fname1, status='Replace', form='unformatted', iostat=info)     
     if (info/=0) then
         DVR0SF = .false.; deallocate(X10, HEig, H2);  return
     end if  
     write(99) N
     write(99) X10(1:N), H2(1:N, 1:N), X1(1:N, 1:N)
     close(99)

     DVR0SF = DVR_STEP2(N, X10, H2, Pot, HEig)     

     if (DVR0SF) then
        open(99, file=fname2, status='Replace', form='unformatted', iostat=info)     
        if (info/=0) then
            DVR0SF = .false.; deallocate(X10, HEig, H2);  return
        end if  
        write(99) N    
        write(99) X10(1:N), HEig(1:N), H2(1:N, 1:N)     
        close(99)

        E0(1:M) = HEIG(1:M)
        DVR0SF = DVR_STEP3(N, M, X10, H2, E0, X0, H0)
     end if  
  end if

  deallocate(X10, HEig, H2) 

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical FUNCTION DVR0RF(POT, M, X0, E0, H0, fname1, fname2)
  implicit none
  integer, intent(in) :: M
  external :: POT
  double precision, intent(out)   :: X0(M), E0(M), H0(M,M)
  character(len=*), intent(IN)    :: fname1, fname2

  logical :: DVR_Step2, DVR_Step3
  double precision, allocatable :: H1(:,:), X10(:), HEig(:)
  integer :: info, N 

  DVR0RF = .false.
   
  open(99, file=fname1, status='Old', form='unformatted', iostat=info)     
  if (info/=0) return   
  read(99) N

  allocate(H1(N,N),X10(N),HEig(N), stat=info)     
  if ((info /= 0) .OR. (N  < M) )then
     close(99); return
  end if

  read(99) X10(1:N), H1(1:N, 1:N)     
  close(99)
     
  DVR0RF = DVR_STEP2(N, X10, H1, Pot, HEig)

  if (DVR0RF) then        
     open(99, file=fname2, status='Replace', form='unformatted', iostat=info)     
     if (info == 0 ) then  
        write(99) N     
        write(99) X10(1:N), HEig(1:N), H1(1:N, 1:N)     
        close(99)

        E0(1:M) = HEIG(1:M)
        DVR0RF = DVR_STEP3(N, M, X10, H1, E0, X0, H0)
     end if
  end if

  deallocate(H1, X10, Heig)

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             Complex         Version                         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical FUNCTION DVRCX0(N, X1, H1, POT, M, X0, E0, H0)
  implicit none
  integer, intent(in) :: N, M
  double complex, intent(inout) :: X1(N,N)
  double precision, intent(IN)  :: H1(N)
  external :: POT
  double precision, intent(out) :: X0(M), E0(M)
  double complex, intent(OUT)   :: H0(M,M)

  logical :: DVRCx_Step0, DVRCx_Step2, DVRCx_Step3
  double precision, allocatable :: X10(:),HEig(:)
  double complex, allocatable   :: H2(:,:)
  integer :: info

  DVRCx0 = .false.;     if ( N < M) return
  allocate(X10(N),HEig(N),H2(N,N),stat=info)
  if (info/=0) return

  DVRCx0 = DVRCx_STEP0(N, X1, H1, X10, H2)

  if (DVRCx0) then

     DVRCx0 = DVRCx_STEP2(N, X10, H2, Pot, HEig)

     if (DVRCx0) then
         E0(1:M) = HEIG(1:M)
         DVRCx0 = DVRCx_STEP3(N, M, X10, H2, E0, X0, H0)
     end if
  end if

  deallocate(X10, HEig, H2) 
end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical FUNCTION DVRCx0SF(N, X1, H1, POT, M, X0, E0, H0, fname1, fname2)
  implicit none
  integer, intent(in) :: N, M
  double precision, intent(inout) :: X1(N,N), H1(N)
  external :: POT
  double precision, intent(out)   :: X0(M), E0(M), H0(M,M)
  character(len=*), intent(IN)    :: fname1, fname2

  logical :: DVRCx_Step0, DVRCx_Step2, DVRCx_Step3
  double precision, allocatable :: X10(:),HEig(:)
  double complex, allocatable   :: H2(:,:)
  integer :: info

  DVRCx0SF = .false.;    if ( N < M) return
  allocate(X10(N),HEig(N),H2(N,N),stat=info)
  if (info/=0) return

  DVRCx0SF = DVRCx_STEP0(N, X1, H1, X10, H2)

  if (DVRCx0SF) then
     open(99, file=fname1, status='Replace', form='unformatted', iostat=info)     
     if (info/=0) then
         DVRCx0SF = .false.; deallocate(X10, HEig, H2); return
     end if
     write(99) N
     write(99) X10(1:N), H2(1:N, 1:N), X1(1:N, 1:N)
     close(99)

     DVRCx0SF = DVRCx_STEP2(N, X10, H2, Pot, HEig)     

     if (DVRCx0SF) then
        open(99, file=fname2, status='Replace', form='unformatted', iostat=info)     
        if (info/=0) then
            DVRCx0SF = .false.; deallocate(X10, HEig, H2); return
        end if  
        write(99) N    
        write(99) X10(1:N), HEig(1:N), H2(1:N, 1:N)     
        close(99)

        E0(1:M) = HEIG(1:M)
        DVRCx0SF = DVRCx_STEP3(N, M, X10, H2, E0, X0, H0)
     end if
  end if

  deallocate(X10, HEig, H2) 

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical FUNCTION DVRCx0RF(POT, M, X0, E0, H0, fname1, fname2)
  implicit none
  integer, intent(in) :: M
  external :: POT
  double precision, intent(out)   :: X0(M), E0(M), H0(M,M)
  character(len=*), intent(IN)    :: fname1, fname2

  logical :: DVRCx_Step2, DVRCx_Step3
  double precision, allocatable ::  X10(:), HEig(:)
  double complex, allocatable :: H1(:,:)
  integer :: info, N 

  DVRCx0RF = .false.
   
  open(99, file=fname1, status='Old', form='unformatted', iostat=info)     
  if (info/=0) return
   
  read(99) N
  allocate(H1(N,N),X10(N),HEig(N), stat=info)     
  if ((info /= 0) .OR. (N  < M) )then
     close(99); return
  end if
  read(99) X10(1:N), H1(1:N, 1:N)     
  close(99)
     
  DVRCx0RF = DVRCx_STEP2(N, X10, H1, Pot, HEig)

  if (DVRCx0RF) then
        
     open(99, file=fname2, status='Replace', form='unformatted', iostat=info)     
     if (info == 0 ) then  
        write(99) N     
        write(99) X10(1:N), HEig(1:N), H1(1:N, 1:N)     
        close(99)

        E0(1:M) = HEIG(1:M)
        DVRCx0RF = DVRCx_STEP3(N, M, X10, H1, E0, X0, H0)
     end if
  end if

  deallocate(H1, X10, Heig)

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
