!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
!c     General subroutine to generate DVR representation          c
!c    logical function DVR_Step1(), DVR_Step2(), DVR_Step3()      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function DVR_Step1(N, X1, H1, X0)
  implicit none
  integer, intent(in) :: N
  double precision, intent(inout) :: X1(N,N), H1(N,N)
  double precision, intent(out)   :: X0(N)

!**************************************************
  include './dvrData.h'
  integer :: lwork, info  
  double precision, allocatable  :: work(:), H1P(:)

!**************************************************
  LWORK = BKSIZE * N ;   dvr_step1 = .false.
  allocate(work(lwork), H1P(N*N),stat=info)

  if (info==0) then 
         ! X1 = V^T*X0*V; X0=diagonal matrix, V=X1
     call  DSYEV('V', 'U', N, X1, N, X0, WORK, LWORK, INFO )   

     if (INFO == 0) then     !H1P=V^T*H1*V
         H1P(1:N*N)=0.0D0
         call dgemm('N','N',N,N,N,1.0D0,H1,N,X1,N,0.0D0,H1P,N)
         call dgemm('T','N',N,N,N,1.0D0,X1,N,H1P,N,0.0D0,H1,N)
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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
!**************************************************
  include './dvrData.h'
  integer :: i, j, symInd, LWORK, INFO  
  double precision    :: work(BKSIZE*N) 
  logical :: symBool
!**************************************************
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  LWORK = BKSIZE * N;  dvr_step2 = .false.
  
  call POT(N, X1, Eig)     ! PES(i) = POT(X0(i))

  do i = 1, N    
     H1(i,i) = H1(i,i) + Eig(i)  
  end do
  
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !c    diagonize H1, get V, and H0 in the second representation      c
  !c     V^T * H0 * V = DIAG_H0,   H0*V = V*DIAG_H0                   c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  call  DSYEV('V', 'U', N, H1, N, Eig, WORK, LWORK, INFO )   ! 99% time spent here
  if (INFO /= 0)   return
 
  ! check all the eigenvectors have the same sign
  symInd = 0;
  do I = 1, N
     if (H1(I,1) /= 0.0D0 ) then
         symBool = .true.
         do J = 1, N
            if (H1(I, J) ==0.0D0) then
                symBool = .false.;  exit
            end if
         end do 
         if (symBool) then
            symInd = i;  exit
         end if
     end if
  end do

  if (symInd==0) return;
  do I=1, N
     if (H1(symInd,I) < 0) H1(1:N, I) = - H1(1:N, I)
  end do

  dvr_step2 = .true.
end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function DVR_Step3(N, M, X1, H1, Eig, X0, H0)
   implicit none
   integer, intent(IN) :: N, M
   double precision, intent(IN) :: X1(N), H1(N,M), Eig(M)
   double precision, intent(OUT):: X0(M), H0(M,M)
 
!**************************************************
   include './dvrData.h'
   integer :: i, j, info, lwork, symInd
   double precision :: X(M,M), work(BKSIZE*M)
   logical :: symBool
!**************************************************

   lwork = BKSIZE*M;  DVR_Step3 = .false.
   if (N<M) return

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !c     calculate X matrix in the second representation      c
  !c                    X = V^T * X0 * V                      c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  do I = 1, M
     do J = I, M
        X(I,J) = sum(H1(1:N,I)*X1(1:N)*H1(1:N,J))
        X(J,I) = X(I, J)
     end do
  end do

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !c  Diagonalize X to get the grid points in PODVR representation  c
  !c        X * V = V *X0,                                          c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  call  DSYEV('V', 'U', M, X, M, X0, WORK, LWORK, INFO )  
  if (INFO /=  0)    return

  ! Check all the eigenvectors have the same sign
  symInd = 0
  do I = 1, M
     if (X(I,1) /= 0.0D0 ) then
         symBool = .true.
         do J = 1, M
            if (X(I, J) ==0.0D0) then
                symBool = .false.;  exit
            end if
         end do
         if (symBool) then
            symInd = i;  exit
         end if
     end if
  end do
  if (symInd==0) return
  do I=1, M
     if (X(symInd,I) < 0) X(1:M, I) = - X(1:M, I)
  end do

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !c        calculate H in the PO DVR representation         c
  !c          H1 = V^T * H0 * V, this is what we need        c
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
  !c    call VTDV(N, DIAG_H0, X, H1)   ! H1 = V^T * H0 * V   c
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  do I = 1, M
     do J = I, M
        H0(I,J) = sum(X(1:M,I)*Eig(1:M)*X(1:M, J))
        H0(J,I) = H0(I, J)
     end do
  end do

  dvr_step3 = .true.

end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

