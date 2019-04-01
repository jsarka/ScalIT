!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
!c     General subroutine to generate DVR representation          c
!c           logical function DVR_Step[0/1/2/3]()                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function DVR_Step0(N, X1, H1, X0, H2)
  implicit none
  integer, intent(in) :: N
  double precision, intent(inout) :: X1(N,N)
  double precision, intent(IN)    :: H1(N)
  double precision, intent(out)   :: X0(N),H2(N,N)

!**************************************************
  include 'dvrData.h'
  integer :: lwork, info, i, j
  double precision, allocatable  :: work(:)
!**************************************************

  lwork = BKSIZE * N ;                dvr_step0 = .false.

  allocate(work(lwork),stat=info)

  if (info==0) then    ! X0 = V^T*X1*V; X0=diagonal matrix

     call  DSYEV('V', 'U', N, X1, N, X0, WORK, lwork, INFO )   

     if (INFO == 0) then     !H2=V^T*H1*V
         do I = 1, N
            do J = 1, N
               H2(I,J) = Sum(X1(1:N,I)*H1(1:N)*X1(1:N,J))
               H2(J,I) = H2(I, J)
            end do
         end do

         DVR_Step0 = .true.
     end if

     deallocate(work)

  end if

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function DVR_Step1(N, X1, H1, X0)
  implicit none
  integer, intent(in) :: N
  double precision, intent(inout) :: X1(N,N), H1(N,N)
  double precision, intent(out)   :: X0(N)

!**************************************************
  include 'dvrData.h'
  integer :: lwork, info
  double precision, allocatable  :: work(:), H1P(:)
!**************************************************

  lwork = BKSIZE * N ;     dvr_step1 = .false.

  allocate(work(lwork), H1P(N*N),stat=info)

  if (info==0) then    ! X0 = V^T*X0*V; X0=diagonal matrix

     call  DSYEV('V', 'U', N, X1, N, X0, work, lwork, info )   

     if (info == 0) then     !H1P=V^T*H1*V
         H1P(1:N*N)=0.0D0
         call DGEMM('N','N',N,N,N,1.0D0,H1,N,X1,N,0.0D0,H1P,N)
         call DGEMM('T','N',N,N,N,1.0D0,X1,N,H1P,N,0.0D0,H1,N)
         DVR_Step1 = .true.
     end if

     deallocate(work,H1P)

  end if

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function DVR_Step2(N, X1, H1, POT, Eig)  
  implicit none
  integer, intent(in) :: N      
  double precision, intent(in)    :: X1(N)
  double precision, intent(inout) :: H1(N, N) 
  external :: POT    
  double precision, intent(out)   :: Eig(N)

!**************************************************
  include 'dvrData.h'
  integer :: i, lwork, info
  double precision, allocatable :: work(:) 
!**************************************************

  lwork = BKSIZE * N;  dvr_step2 = .false.
  
  allocate(work(lwork),stat=info)

  if (info==0) then
      call POT(N, X1, Eig)     ! PES(i) = POT(X0(i))

      do i = 1, N    
         H1(i,i) = H1(i,i) + Eig(i)  
      end do  

      call  DSYEV('V', 'U', N, H1, N, Eig, WORK, lwork, info )
 
      dvr_step2 = (info==0)
     
      deallocate(work)

   end if

end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function DVR_Step3(N, M, X1, H1, Eig, X0, H0)
   implicit none
   integer, intent(IN) :: N, M
   double precision, intent(IN) :: X1(N), H1(N,M), Eig(M)
   double precision, intent(OUT):: X0(M), H0(M,M)
 
!**************************************************
   include 'dvrData.h'
   integer :: i, j,  info, lwork
   double precision :: X(M,M), work(BKSIZE*M)
!**************************************************

   lwork = BKSIZE*M;  DVR_Step3 = .false.

   if (N < M) return

           ! X = V^T*X1*V
   do I = 1, M
      do J = I, M
         X(I,J) = sum(H1(1:N,I)*X1(1:N)*H1(1:N,J))
         X(J,I) = X(I, J)
      end do
   end do

   call  DSYEV('V', 'U', M, X, M, X0, WORK, lwork, info )  

   if (info ==  0)  then   ! H0 = V^T * Eig * V 
      do I = 1, M         
         do J = I, M
            H0(I,J) = sum(X(1:M,I)*Eig(1:M)*X(1:M, J))
            H0(J,I) = H0(I, J)
         end do
      end do

      dvr_step3 = .true.
  end if

end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

