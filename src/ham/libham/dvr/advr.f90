!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      Non-orthogonal basis set for DVR                  c
!c  S1=overlap matrx <f|f>, X1=x matrix <f|x|f>           c
!c  S1=V*S0^(-1/2)*U, X1=S1^-1=U^H*S0^(1/2)*V^H           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function ADVR(N, S1, X1, X0)
   implicit none
   include 'dvrData.h'
   integer, intent(IN) :: N
   double precision, intent(INOUT) :: S1(N,N), X1(N,N)
   double precision, intent(OUT) :: X0(N)

   double precision, allocatable :: X2(:,:), tmp(:,:), work(:), Y0(:)
   integer :: info, i, lwork
   
   ADVR = .false.;     lwork=BKSIZE*N
   allocate(X2(N,N), tmp(N,N), work(lwork), Y0(N), stat=info)   
   if (info/=0)  return

   call  DSYEV('V', 'U', N, S1, N, Y0, work, lwork, info) 
  
   if ((info/=0).or. (Y0(1)<=0.0D0)) then
       deallocate(X2, tmp, work, Y0);  return
   end if

   Y0(1:N)=sqrt(Y0(1:N))

   do i = 1, N
      tmp(1:N, i) = S1(1:N,i)/Y0(i)
   end do
   call VTAV(N, N, X1, tmp, X2)   
   call  DSYEV('V', 'U', N, X2, N, X0, work, lwork, info) 

   if (info==0) then
      do i = 1, N
          tmp(1:N, i) = S1(1:N,i)*Y0(i)
      end do

      call DGEMM('T','T', N, N, N,1.0D0,X2,N,tmp,N,0.0D0,X1,N)

      do i = 1, N
          tmp(1:N, i) = S1(1:N,i)/Y0(i)
      end do

      call DGEMM('N','N', N, N, N,1.0D0,tmp,N,X2,N,0.0D0,S1,N)

      ADVR = .true.
   end if

   deallocate(X2, tmp, work, Y0);

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      Non-orthogonal basis set for DVR                  c
!c  S1=overlap matrx <f|f>, X1=x matrix <f|x|f>           c
!c  S1=V*S0^(-1/2)*U, X1=S1^-1=U^H*S0^(1/2)*V^H           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function ADVRF(N, S1, X1, X0, wrFlag, fname)
   implicit none
   include 'dvrData.h'
   integer, intent(IN) :: N
   double precision, intent(INOUT) :: S1(N,N), X1(N,N)
   double precision, intent(OUT) :: X0(N)
   logical, intent(IN) :: wrFlag
   character(len=*), intent(IN) :: fname

   double precision, allocatable :: X2(:,:), tmp(:,:), work(:), Y0(:)
   integer :: info, i, lwork, N0
   
   ADVRF = .false.;     lwork=BKSIZE*N
   allocate(X2(N,N), tmp(N,N), work(lwork), Y0(N), stat=info)   
   if (info/=0)  return

   if (wrFlag) then    ! save transformation matrices in file
      call  DSYEV('V', 'U', N, S1, N, Y0, work, lwork, info) 
  
      if ((info/=0).or. (Y0(1)<=0.0D0)) then
          deallocate(X2, tmp, work, Y0);  return
      end if
      Y0(1:N)=sqrt(Y0(1:N))

      do i = 1, N
         tmp(1:N, i) = S1(1:N,i)/Y0(i)
      end do
      call VTAV(N, N, X1, tmp, X2)
      call  DSYEV('V', 'U', N, X2, N, X0, work, lwork, info) 

      if (info==0) then

          open(99,file=fname,status='Replace',form='unformatted',iostat=info)

          if (info==0) then
              write(99) N
              write(99) Y0(1:N), X0(1:N)
              write(99) S1(1:N,1:N), X2(1:N,1:N)                
          end if

          close(99);  ADVRF=(info==0)

      end if
   else  ! load data from file
      open(99,file=fname,status='Old',form='unformatted',iostat=info)

      if (info==0) then
          read(99) N0

          if (N0==N) then
              read(99) Y0(1:N), X0(1:N)
              read(99) S1(1:N,1:N), X2(1:N,1:N)                          
          end if       

          close(99);       ADVRF = (N0==N)
      end if
   end if 

   if (ADVRF) then
      do i = 1, N
         tmp(1:N, i) = S1(1:N,i)*Y0(i)
      end do

      call DGEMM('T','T', N, N, N,1.0D0,X2,N,tmp,N,0.0D0,X1,N)

      do i = 1, N
         tmp(1:N, i) = S1(1:N,i)/Y0(i)
      end do

      call DGEMM('N','N', N, N, N,1.0D0,tmp,N,X2,N,0.0D0,S1,N)      
   end if

   deallocate(X2, tmp, work, Y0);

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function ADVRCx(N, S1Cx, X1Cx, X0)
   implicit none
   include 'dvrData.h'
   integer, intent(IN) :: N
   double complex, intent(INOUT) :: S1Cx(N,N), X1Cx(N,N)
   double precision, intent(OUT) :: X0(N)

   double complex,parameter :: alpha=(1.0D0,0.0D0), beta=(0.0D0,0.0D0)
   double complex,   allocatable :: X2Cx(:,:), tmp(:,:), work(:),rwork(:)
   double precision, allocatable :: Y0(:)
   integer :: info, i, lwork
   
   ADVRCx = .false.;     lwork=BKSIZE*N
   allocate(X2Cx(N,N),tmp(N,N),work(lwork),rwork(lwork),Y0(N),stat=info)   
   if (info/=0)  return

   call  ZHEEV('V', 'U', N, S1Cx, N, Y0, work, lwork, rwork, info) 
  
   if ((info/=0) .OR. (Y0(1)<=0.0D0)) then
       deallocate(X2Cx, tmp, work, rwork, Y0);  return
   end if
   Y0(1:N)=sqrt(Y0(1:N))

   do i = 1, N
      tmp(1:N, i) = S1Cx(1:N,i)/Y0(i)
   end do

   call VHAVCx(N, N, X1Cx, tmp, X2Cx)   
   call  ZHEEV('V', 'U', N, X2Cx, N, X0, work, lwork, rwork, info) 

   if (info==0) then
      do i = 1, N
          tmp(1:N, i) = S1Cx(1:N,i)*Y0(i)
      end do

      call ZGEMM('C','C', N, N, N,alpha,X2Cx,N,tmp,N,beta,X1Cx,N)

      do i = 1, N
          tmp(1:N, i) = S1Cx(1:N,i)/Y0(i)
      end do

      call ZGEMM('N','N', N, N, N,alpha,tmp,N,X2Cx,N,beta,S1Cx,N)

      ADVRCx = .true.
   end if

   deallocate(X2Cx, tmp, work, rwork, Y0);

end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function ADVRCxF(N, S1Cx, X1Cx, X0, wrFlag, fname)
   implicit none
   include 'dvrData.h'
   integer, intent(IN) :: N
   double complex, intent(INOUT) :: S1Cx(N,N), X1Cx(N,N)
   double precision, intent(OUT) :: X0(N)
   logical, intent(IN) :: wrFlag
   character(len=*), intent(IN) :: fname

   double complex,parameter :: alpha=(1.0D0,0.0D0), beta=(0.0D0,0.0D0)
   double complex,   allocatable :: X2Cx(:,:), tmp(:,:), work(:),rwork(:)
   double precision, allocatable :: Y0(:)
   integer :: info, i, lwork, N0
   
   ADVRCxF = .false.;     lwork=BKSIZE*N;  
   allocate(X2Cx(N,N),tmp(N,N),work(lwork),rwork(lwork),Y0(N),stat=info)   
   if (info/=0)  return

   if (wrFlag) then    ! save transformation matrices in fil
       call  ZHEEV('V', 'U', N, S1Cx, N, Y0, work, lwork, rwork, info) 
  
       if ((info/=0).or. (Y0(1)<=0.0D0)) then
           deallocate(X2Cx, tmp, work, rwork, Y0);  return
       end if

       Y0(1:N)=sqrt(Y0(1:N))

       do i = 1, N
          tmp(1:N, i) = S1Cx(1:N,i)/Y0(i)
       end do
       call VHAVCx(N, N, X1Cx, tmp, X2Cx)   
       call  ZHEEV('V', 'U', N, X2Cx, N, X0, work, lwork, rwork, info) 

       if (info==0) then          ! save matrices
          open(99,file=fname,status='Replace',form='unformatted',iostat=info)

          if (info==0) then
              write(99) N
              write(99) Y0(1:N), X0(1:N)
              write(99) S1Cx(1:N,1:N), X2Cx(1:N,1:N)              
          end if

          close(99);   ADVRCxF=(info==0)
       end if
   else
      open(99,file=fname,status='Old',form='unformatted',iostat=info)

      if (info==0) then
          read(99) N0

          if (N0==N) then
              read(99) Y0(1:N), X0(1:N)
              read(99) S1Cx(1:N,1:N), X2Cx(1:N,1:N)                          
          end if       

          close(99);    ADVRCxF=(N==N0)
       end if
   end if

   if (ADVRCxF) then
      do i = 1, N
         tmp(1:N, i) = S1Cx(1:N,i)*Y0(i)
      end do

      call ZGEMM('C','C', N, N, N,alpha,X2Cx,N,tmp,N,beta,X1Cx,N)

      do i = 1, N
         tmp(1:N, i) = S1Cx(1:N,i)/Y0(i)
      end do

      call ZGEMM('T','T', N, N, N,alpha,tmp,N,X2Cx,N,beta,S1Cx,N)
   end if

   deallocate(X2Cx, tmp, work, rwork, Y0);
end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
