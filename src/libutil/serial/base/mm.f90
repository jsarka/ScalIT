!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                B = V^T*H*V                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine VTHV(N, M, H, V, B)   
    implicit none
    integer, intent(IN) :: N, M
    double precision, intent(IN)  :: H(N,N), V(N,M)
    double precision, intent(OUT) :: B(M,M)

    double precision :: tmp(N,M)
 
    call dgemm('N','N',N,M,N,1.0D0,H,N,V,N,0.0D0,tmp,N)

    call dgemm('T','N',M,M,N,1.0D0,V,N,tmp,N,0.0D0,B,M)

  end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     H = V^T * D * V,  D is a diagonal matrix      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine VTDV(N, M, D, V, H)
  implicit none
  integer, intent(IN) :: N,M
  double precision, intent(IN) :: D(N), V(N,M)
  double precision, intent(OUT) :: H(M,M) 

  integer :: i, j

  do i = 1, M
     do j = 1, I
        H(I, J) = SUM(V(1:N,I)*D(1:N)*V(1:N, J))
        H(J, I) = H(I, J)
     end do
  end do        

end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc



!************************************************
!*        Diagonalize a Symmetric Matrix        *
!************************************************
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
LOGICAL FUNCTION DIAG(JOBZ, N, H1, E0)
   IMPLICIT NONE
   CHARACTER,INTENT(IN):: JOBZ
   INTEGER, INTENT(IN) :: N
   DOUBLE PRECISION, INTENT(INOUT) :: H1(N,N)
   DOUBLE PRECISION, INTENT(OUT)   :: E0(N)

   INTEGER :: INFO, lwork, symInd, I, J
   logical :: symBool
   integer, parameter :: scale=3
   double precision   :: WORK(scale*N) 
   
   lwork = scale*N;    diag = .false.
   CALL DSYEV(JOBZ,'U', N, H1, N, E0, WORK, LWORK, Info)

   if (Info /= 0) return

   ! Make sure all the eigenvectors have the same sign
   symInd = 0;
   DO I = 1, N
     if (H1(I,1) /= 0.0D0 ) then
         symBool = .true.
         DO J = 1, N
            IF (H1(I, J) ==0.0D0) THEN
                symBool = .false.;  EXIT
            END IF
         END DO 
         IF (symBool) THEN
            symInd = i;  EXIT
         END IF
     end if
  END DO
  IF (symInd==0) return;
  DO I=1, N
     if (H1(symInd,I) < 0) H1(1:N, I) = - H1(1:N, I)
  END DO

  diag = .TRUE.

END FUNCTION
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c          Compact the matrix                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function CompactMat(E0, N, H1, M, H2)
   implicit none
   double precision, intent(IN)    :: E0
   integer, intent(IN)             :: N
   double precision, intent(INOUT) :: H1(N,N)
   integer, intent(OUT)            :: M
   double precision, intent(OUT)   :: H2(N,N)
   
   integer, parameter :: SC=3
   integer :: lwork, info, i
   double precision, allocatable :: HT(:,:), work(:), H0(:)
   double precision, allocatable :: t1(:,:),t2(:,:)

   lwork = SC * N
   CompactMat = .false.
   allocate(HT(N,N), work(lwork), H0(N), stat=info) 
   if (info/=0) return  

   H2(1:N,1:N) = H1(1:N,1:N)

   ! Diagonalize H1
   call DSYEV('V','U',N,H1,N,H0,work,lwork,info)
   if (info/=0) then
       deallocate(HT, work, H0);   return
   end if
   
   print *, 'Eigenvalue:', H0

   compactMat = .true.
   if (H0(N)<=E0) then  ! No compacting
       M = N;  deallocate(HT, work, H0); return 
   end if
   if (H0(1)>E0) then
       M = 0;  deallocate(HT, work, H0); return
   end if

   ! compacting, H2(M,M) = V^T(M,N)*E0(N)I*V(N,M)
   !  V=H1, H1=H2
   M = 0
   do i = 1, N
      if (H0(i)>E0) then
          M=i-1; exit
      end if
   end do

   allocate(t1(N,M),t2(M,M))
!   t1 = MatMul(H2(1:N,1:N),H1(1:N,1:M))
   do i=1,N
      t1(i,1:M) = H0(i)*H1(i,1:M)
   end do
   H1 = transpose(H1)
   t2 = MatMul(H1(1:M,1:N),t1(1:N,1:M))   
   H2(1:M,1:M)=t2(1:M,1:M)
   deallocate(t1,t2)  

!   call DGEMM('N','N',N,M,N,1.0D0,H2,N,H1,N,0.0D0,HT,N)
!   call DGEMM('T','N',M,M,N,1.0D0,H1,N,HT,N,0.0D0,H2,N)

   deallocate(HT, work, H0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
