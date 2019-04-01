!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         Fill the Diagonal elements of matrix      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine fillDiagMat_DX(NIN, NOUT, diagMat, Mat)
     implicit none
     integer, intent(IN) :: NIN, NOUT
     double precision, intent(IN) :: diagMat(NIN, NOUT, NOUT)
     double complex, intent(OUT)  :: mat(NIN*NOUT, NIN*NOUT)

     integer :: i, j, k
     integer :: rowInd, colInd

     mat(1:NIN*NOUT, 1:NIN*NOUT) = 0.0D0
     do I = 1, NOUT         ! column
        do J = 1, NIN
           colInd = (I-1)*NIN + J
           do K = 1, NOUT
               rowInd = (k-1)*NIN + J
               mat(rowInd, colInd) = diagMat(J, K, I)
           end do
        end do
     end do

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine addDiagMat_DX(NIN, NOUT, diagMat, Mat)
     implicit none
     integer, intent(IN) :: NIN, NOUT
     double precision, intent(IN)  :: diagMat(NIN, NOUT, NOUT) 
     double complex, intent(INOUT) :: mat(NIN*NOUT, NIN*NOUT)

     integer :: i, j, k
     integer :: rowInd, colInd

     do I = 1, NOUT         ! column
        do J = 1, NIN
           colInd = (I-1)*NIN + J
           do K = 1, NOUT
               rowInd = (k-1)*NIN + J
               mat(rowInd, colInd) = mat(rowInd, colInd)+diagMat(J, K, I)
           end do
        end do
     end do

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Matrix-vector product used in OSB package:             c
!c                     XOUT = A * X                           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine H1X_DX(NIN, NOUT, H1, X, XOUT)
  implicit none  
  integer,intent(IN) :: NIN, NOUT                
  double precision, intent(IN)  :: H1(NOUT,NOUT)
  double complex, intent(IN)   :: X(NIN,NOUT)
  double complex, intent(OUT)  :: XOUT(NIN,NOUT)

  integer :: i,j

  do i = 1, NOUT
     XOUT(1:NIN, i) = H1(i, 1)*X(1:NIN, 1) 
     do j = 2, NOUT
         XOUT(1:NIN, i) = XOUT(1:NIN, i) + H1(i, j)*X(1:NIN, j) 
     end do
  end do
end subroutine 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Matrix-vector product used in OSB package:             c
!c                     XOUT = A^T * X                         c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine  H1TX_DX(NIN, NOUT, H1, X, XOUT)
  implicit none
  integer,intent(IN) :: NIN, NOUT                
  double precision,intent(IN) :: H1(NOUT,NOUT)
  double complex,intent(IN) :: X(NIN,NOUT)
  double complex,intent(OUT):: XOUT(NIN,NOUT)

  integer ::    I,J

  do i = 1, NOUT
     XOUT(1:NIN, i) = H1(1, i)*X(1:NIN, 1) 
     do j = 2, NOUT
         XOUT(1:NIN, i) = XOUT(1:NIN, i) + H1(j, i)*X(1:NIN, j) 
     end do
  end do
end subroutine 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Matrix-vector product used in OSB package:             c
!c                     XOUT = A * X                           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine H2X_DX(NIN, NOUT, H1, CORX, X, XOUT)
  implicit none 
  integer,intent(IN) :: NIN, NOUT                
  double precision, intent(IN) :: H1(NOUT,NOUT)
  double precision, intent(IN) :: CORX(NIN)
  double complex, intent(IN)   :: X(NIN,NOUT)
  double complex, intent(OUT)  :: XOUT(NIN,NOUT)

  integer :: i,j,k

  do i=1,NOUT
     XOUT(1:NIN, i) = H1(i,1) * CORX(1:NIN) * X(1:NIN, 1)
     do j=2,NOUT
         XOUT(1:NIN,i) = XOUT(1:NIN,i)+H1(i,j)*CORX(1:NIN)*X(1:NIN,j)
     end do
  end do
end subroutine 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Matrix-vector product used in OSB package:             c
!c                     XOUT = A * X                           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine H2TX_DX(NIN, NOUT, H1, CORX, X, XOUT)
  implicit none 
  integer,intent(IN) :: NIN, NOUT                
  double precision, intent(IN) :: H1(NOUT,NOUT)
  double precision, intent(IN) :: CORX(NIN)
  double complex, intent(IN)   :: X(NIN,NOUT)
  double complex, intent(OUT)  :: XOUT(NIN,NOUT)

  integer :: i,j

  do i=1,NOUT
     XOUT(1:NIN, i) = H1(1, i) * CORX(1:NIN) * X(1:NIN, 1)
     do j=2,NOUT
         XOUT(1:NIN,i) = XOUT(1:NIN,i)+H1(j,i)*CORX(1:NIN)*X(1:NIN,j)
     end do
  end do
end subroutine 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Matrix-vector product used in OSB package:             c
!c                     XOUT = A * X                           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine H3X_DX(NIN, NOUT, A, X, XOUT) 
  integer,intent(IN) ::  NIN,NOUT
  double precision, intent(IN) :: A(NIN,NOUT,NOUT)
  double complex, intent(IN)   :: X(NIN, NOUT)
  double complex, intent(OUT)  :: XOUT(NIN, NOUT)

  integer  ::  i, j
      
  do i=1,NOUT     
     XOUT(1:NIN, i) = A(1:NIN, i, 1) * X(1:NIN, 1)
     do j = 2, NOUT
           XOUT(1:NIN, i) = XOUT(1:NIN, i) + A(1:NIN, i, j) * X(1:NIN, j)
     end do   
  end do
  
end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Matrix-vector product used in OSB package:             c
!c                     XOUT = A * X                           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine H3TX_DX(NIN, NOUT, A, X, XOUT) 
  implicit none
  integer,intent(IN) ::  NIN,NOUT
  double precision, intent(IN) :: A(NIN, NOUT,NOUT)
  double complex, intent(IN)   :: X(NIN, NOUT)
  double complex, intent(OUT)  :: XOUT(NIN, NOUT)

  integer  ::  i, j
      
  do i=1, NOUT     
     XOUT(1:NIN, i) = A(1:NIN, 1, i) * X(1:NIN, 1)
     do j = 2, NOUT
           XOUT(1:NIN, i) = XOUT(1:NIN, i) + A(1:NIN, j, i) * X(1:NIN, j)
     end do   
  end do  
end subroutine 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Calculate H = V^T * D * V,  D is diagonal matrix     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   H(i, j) = sum(k) [V(k,i) * d(k) * V(k, j)]           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine VTDV_DX(N, D, V, H)
  implicit none
  integer, intent(in)   :: N
  double complex, intent(in)   :: D(N)
  double precision, intent(in) :: V(N,N)
  double complex, intent(out)  :: H(N,N)

  integer :: i, j, k
  double complex :: tmp

  do i = 1, N
     do j = 1, N
        H(i, j) = sum( V(1:N,i) * D(1:N) * V(1:N,j) )
     end do
  end do        
end subroutine 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Calculate H1 = V^T * H0 * V,  D is a normal matrix    c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   H(i, j) = sum(k) [V(k,i) * H0(k, l) * V(k, j)]       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine VTHV_DX(N, H0, V, H1)
  implicit none
  integer, intent(in)   :: N
  double complex, intent(in)   :: H0(N,N)
  double precision, intent(in) :: V(N,N)
  double complex, intent(out)  :: H1(N,N)

  integer :: i, j, k
  double complex :: tmp

  do i = 1, N
     do j = 1, N
        tmp = 0.0D0
        do k = 1, N
           tmp = tmp + sum(V(1:N, i) * H0(1:N, k)) * V(k, j)
        end do
        H1(i, j) = tmp
     end do
  end do        
 
end subroutine 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
