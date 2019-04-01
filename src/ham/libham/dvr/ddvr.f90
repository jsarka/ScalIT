!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Prepare non-orthogonal basis sets for DVR         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         S0=V^T*S1*V; S0^(1/2)=>S0, V=>S1              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function DDVR(N, S1, S0)    
   implicit none
   include 'dvrData.h'
   integer, intent(IN) :: N
   double precision, intent(INOUT) :: S1(N,N)
   double precision, intent(OUT)   :: S0(N)

!**************************************************
   integer :: lwork, info
   double precision, allocatable  :: work(:)
!**************************************************

   lwork = BKSIZE * N ;                DDVR = .false.

   allocate(work(lwork),stat=info)

   if (info==0) then  

      call  DSYEV('V', 'U', N, S1, N, S0, work, lwork, info)
  
      deallocate(work)

      DDVR = (info==0)

  end if

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         S0=V^T*S1*V; S0^(1/2)=>S0, V=>S1              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function DDVRCX(N, S1Cx, S0)    
   implicit none
   include 'dvrData.h'
   integer, intent(IN) :: N
   double complex, intent(INOUT) :: S1Cx(N,N)
   double precision, intent(OUT) :: S0(N)

!**************************************************
   integer :: lwork, info
   double complex, allocatable  :: work(:), rwork(:)
!**************************************************

   lwork = BKSIZE * N ;            DDVRCX = .false.

   allocate(work(lwork),rwork(lwork),stat=info)

   if (info==0) then  

      call  ZHEEV('V', 'U', N, S1Cx, N, S0, work, lwork, rwork, info)
  
      deallocate(work, rwork)

      DDVRCX = (info==0)

   end if

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         S0=V^T*S1*V; S0^(1/2)=>S0, V=>S1              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function DDVRF(N, S1, S0, wrFlag, fname)    
   implicit none
   integer, intent(IN) :: N
   double precision, intent(INOUT) :: S1(N,N)
   double precision, intent(OUT)   :: S0(N)
   logical, intent(IN) :: wrFlag
   character(len=*), intent(IN) :: fname

!**************************************************
   integer :: n0, info
   logical :: DDVR
!**************************************************

   DDVRF = .false.

   if (wrFlag) then    ! save data

      if (DDVR(N,S1,S0)) then 

         open(99, file=fname, status='Replace', form='unformatted', iostat=info)
     
         if (info/=0) return     
  
         write(99) N

         write(99) S0(1:N), S1(1:N, 1:N)
 
         close(99)
        
         DDVRF = .TRUE.

      end if

   else     ! read data      

      open(99, file=fname, status='Old', form='unformatted', iostat=info)
     
      if (info/=0) return     
  
      read(99) N0
      
      if (N0==N)    read(99) S0(1:N), S1(1:N, 1:N)

      close(99)

      DDVRF = (N0==N)

   end if

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         S0=V^T*S1*V; S0^(1/2)=>S0, V=>S1              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function DDVRCXF(N, S1Cx, S0, wrFlag, fname)    
   implicit none
   integer, intent(IN) :: N
   double complex, intent(INOUT) :: S1Cx(N,N)
   double precision, intent(OUT) :: S0(N)
   logical, intent(IN) :: wrFlag
   character(len=*), intent(IN) :: fname

!**************************************************
   integer :: n0, info
   logical :: DDVRCX
!**************************************************

   DDVRCXF = .false.

   if (wrFlag) then    ! save data

      if (DDVRCX(N,S1Cx,S0)) then 

         open(99, file=fname, status='Replace', form='unformatted', iostat=info)
     
         if (info/=0) return     
  
         write(99) N

         write(99) S0(1:N), S1Cx(1:N, 1:N)
 
         close(99)
        
         DDVRCXF = .TRUE.

      end if

   else     ! read data      

      open(99, file=fname, status='Old', form='unformatted', iostat=info)
     
      if (info/=0) return     
  
      read(99) N0
      
      if (N0==N)    read(99) S0(1:N), S1Cx(1:N, 1:N)

      close(99)

      DDVRCXF = (N0==N)

   end if

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
