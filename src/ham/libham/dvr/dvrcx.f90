!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
!c     General subroutine to generate DVR representation          c
!c    logical function DVR_Step1(), DVR_Step2(), DVR_Step3()      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function DVRCx_Step0(N, X1, H1, X0, H2)
  implicit none
  integer, intent(in) :: N
  double complex, intent(inout) :: X1(N,N)
  double precision, intent(IN)  :: H1(N)
  double precision, intent(out) :: X0(N)
  double complex, intent(OUT)   :: H2(N,N)

!**************************************************
  include 'dvrData.h'  
  integer :: lwork, info, i, j
  double complex, allocatable  :: rwork(:), work(:)

!**************************************************
  lwork = BKSIZE * N ;        DVRCx_step0 = .false.

  allocate(rwork(lwork),work(lwork),stat=info)

  if (info==0) then 

     call  ZHEEV('V', 'U', N, X1, N, X0, WORK, LWORK, RWork, INFO )   

     if (info == 0) then     !H2=V^T*H1*V

         do I = 1, N
            do J = 1, N
               H2(I,J) = Sum(Conjg(X1(1:N,I))*H1(1:N)*X1(1:N,J))
               H2(J,I) = Conjg(H2(I, J))
            end do
         end do

         DVRCx_Step0 = .true.

     end if

     deallocate(work, rwork)

  end if

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function DVRCx_Step1(N, X1, H1, X0)
  implicit none
  integer, intent(in) :: N
  double complex, intent(inout) :: X1(N,N), H1(N,N)
  double precision, intent(out) :: X0(N)

!**************************************************
  include 'dvrData.h'
  double complex, parameter :: ONECx=(1.0D0,0.0D0),ZEROCX=(0.0D0,0.0D0)
  integer :: lwork, info
  double complex, allocatable  :: work(:), H1P(:), rwork(:)

!**************************************************
  lwork = BKSIZE * N ;   DVRCx_step1 = .false.

  allocate(work(lwork), rwork(lwork), H1P(N*N), stat=info)

  if (info==0) then 

     call  ZHEEV('V', 'U', N, X1, N, X0, WORK, lwork, rwork, info )   

     if (info == 0) then     !H1P=V^T*H1*V
         H1P(1:N*N)=0.0D0
         call ZGemm('N','N',N,N,N,ONECX,H1,N,X1,N,ZEROCX,H1P,N)
         call ZGemm('C','N',N,N,N,ONECX,X1,N,H1P,N,ZEROCX,H1,N)

         DVRCx_Step1 = .true.

     end if

     deallocate(work,rwork,H1P)

  end if

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function DVRCx_Step2(N, X1, H1, POT, Eig)  
  implicit none
  integer, intent(in) :: N      
  double precision, intent(in)  :: X1(N)
  double complex, intent(inout) :: H1(N, N) 
  external :: POT    
  double precision, intent(out) :: Eig(N)

!**************************************************
  include 'dvrData.h'
  integer :: i, lwork, info
  double complex, allocatable :: work(:), rwork(:)
!**************************************************

  lwork = BKSIZE * N;  DVRCx_step2 = .false.
  allocate(work(lwork), rwork(lwork),stat=info)

  if (info==0) then  

      call POT(N, X1, Eig)     ! PES(i) = POT(X0(i))

      do i = 1, N    
         H1(i,i) = H1(i,i) + Eig(i)  
      end do

      call  ZHEEV('V', 'U', N, H1, N, Eig, work, lwork, rwork, info )   
  
      DVRCx_step2 = (info==0)
    
      deallocate(work, rwork)
  
  end if

!  print *, 'Eig1:',Eig
!  print *, 'H1:',H1

end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function DVRCx_Step3(N, M, X1, H1, Eig, X0, H0)
   implicit none
   integer, intent(IN) :: N, M
   double precision, intent(IN) :: X1(N), Eig(M)
   double complex, intent(IN)   :: H1(N,M)
   double precision, intent(OUT):: X0(M)
   double complex, intent(OUT)  :: H0(M,M)
 
!**************************************************
   include 'dvrData.h'
   integer :: i, j,  info, lwork
   double complex :: X(M,M), work(BKSIZE*M), rwork(BKSIZE*M)
!**************************************************

   lwork = BKSIZE*M;  DVRCx_Step3 = .false.

   if ( N < M ) return

   do I = 1, M           ! X = V^T * X0 * V
      do J = I, M
         X(I,J) = Sum(Conjg(H1(1:N,I))*X1(1:N)*H1(1:N,J))
         X(J,I) = Conjg(X(I, J))
      end do
   end do

   call  ZHEEV('V', 'U', M, X, M, X0, WORK, lwork, rwork, info )  

   if (info ==  0) then   ! H0=V^H*Eig*V

      do I = 1, M
         do J = I, M
            H0(I,J) = Sum(Conjg(X(1:M,I))*Eig(1:M)*X(1:M, J))
            H0(J,I) = Conjg(H0(I, J))
         end do
      end do

      DVRCx_step3 = .true.

   end if

end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

