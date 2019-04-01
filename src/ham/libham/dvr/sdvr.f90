!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Prepare non-orthogonal basis sets for DVR         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         S0=V^T*S1*V; S0^(1/2)=>S0, V=>S1              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function SDVR(N, S1, S0)    
   implicit none
   include 'dvrData.h'
   integer, intent(IN) :: N
   double precision, intent(INOUT) :: S1(N,N)
   double precision, intent(OUT)   :: S0(N)

!**************************************************
   integer :: lwork, info
   double precision, allocatable  :: work(:)
!**************************************************

   lwork = BKSIZE * N ;                SDVR = .false.

   allocate(work(lwork),stat=info)

   if (info==0) then  

      call  DSYEV('V', 'U', N, S1, N, S0, work, lwork, info)
  
      deallocate(work)

      SDVR = (info==0)

      if (SDVR)  S0(1:N) = DSQRT(S0(1:N))

  end if

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         S0=V^T*S1*V; S0^(1/2)=>S0, V=>S1              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function SDVRCX(N, S1, S0)    
   implicit none
   include 'dvrData.h'
   integer, intent(IN) :: N
   double complex, intent(INOUT) :: S1(N,N)
   double precision, intent(OUT) :: S0(N)

!**************************************************
   integer :: lwork, info
   double complex, allocatable  :: work(:), rwork(:)
!**************************************************

   lwork = BKSIZE * N ;                SDVRCX = .false.

   allocate(work(lwork),rwork(lwork),stat=info)

   if (info==0) then  

      call  ZHEEV('V', 'U', N, S1, N, S0, work, lwork, rwork, info)
  
      deallocate(work, rwork)

      SDVRCX = (info==0)

      if (SDVRCX)  S0(1:N) = DSQRT(S0(1:N))

   end if

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         S0=V^T*S1*V; S0^(1/2)=>S0, V=>S1              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function SDVRF(N, S1, S0, wrFlag, fname)    
   implicit none
   integer, intent(IN) :: N
   double precision, intent(INOUT) :: S1(N,N)
   double precision, intent(OUT)   :: S0(N)
   logical, intent(IN) :: wrFlag
   character(len=*), intent(IN) :: fname

!**************************************************
   integer :: n0, info
   logical :: SDVR
!**************************************************

   SDVRF = .false.

   if (wrFlag) then    ! save data

      if (SDVR(N,S1,S0)) then 

         open(99, file=fname, status='Replace', form='unformatted', iostat=info)
     
         if (info/=0) return     
  
         write(99) N

         write(99) S0(1:N), S1(1:N, 1:N)
 
         close(99)
        
         SDVRF = .TRUE.

      end if

   else     ! read data      

      open(99, file=fname, status='Old', form='unformatted', iostat=info)
     
      if (info/=0) return     
  
      read(99) N0
      
      if (N0==N)    read(99) S0(1:N), S1(1:N, 1:N)

      close(99)

      SDVRF = (N0==N)

   end if

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         S0=V^T*S1*V; S0^(1/2)=>S0, V=>S1              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function SDVRCXF(N, S1, S0, wrFlag, fname)    
   implicit none
   integer, intent(IN) :: N
   double complex, intent(INOUT) :: S1(N,N)
   double precision, intent(OUT)   :: S0(N)
   logical, intent(IN) :: wrFlag
   character(len=*), intent(IN) :: fname

!**************************************************
   integer :: n0, info
   logical :: SDVRCX
!**************************************************

   SDVRCXF = .false.

   if (wrFlag) then    ! save data

      if (SDVRCX(N,S1,S0)) then 

         open(99, file=fname, status='Replace', form='unformatted', iostat=info)
     
         if (info/=0) return     
  
         write(99) N

         write(99) S0(1:N), S1(1:N, 1:N)
 
         close(99)
        
         SDVRCXF = .TRUE.

      end if

   else     ! read data      

      open(99, file=fname, status='Old', form='unformatted', iostat=info)
     
      if (info/=0) return     
  
      read(99) N0
      
      if (N0==N)    read(99) S0(1:N), S1(1:N, 1:N)

      close(99)

      SDVRCXF =  (N0==N)

   end if

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
