!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                    makeI                               c
!c   Create a Identity matrix: A(i,i)=1, A(i,j) = 0       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine makeI_CX(N, A)
  implicit none
  integer, intent(IN) :: N
  double complex, intent(OUT) :: A(N, N)

  integer :: i

  A(1:N, 1:N) = 0.0D0
  do  i = 1, N
     A(i, i)   = 1.0D0
  end do
end subroutine 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                Randomly generated matrix                c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine randMat_CX(n,m,mat)
    implicit none
    integer, intent(IN) :: N, M
    double complex, intent(OUT) :: mat(N,M)    

    double precision, dimension(N, M) :: A, B

    call random_number(A)
    call random_number(B)

    mat(1:N, 1:M) = DCMPLX(A(1:N, 1:M), B(1:N, 1:M))
end subroutine 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         Fill the Diagonal elements of matrix      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine fillDiagMat_CX(NIN, NOUT, diagMat, Mat)
     implicit none
     integer, intent(IN) :: NIN, NOUT
     double complex, intent(IN) :: diagMat(NIN, NOUT, NOUT)
     double complex, intent(OUT):: mat(NIN*NOUT, NIN*NOUT)

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
subroutine addDiagMat_CX(NIN, NOUT, diagMat, Mat)
     implicit none
     integer, intent(IN) :: NIN, NOUT
     double complex, intent(IN) :: diagMat(NIN, NOUT, NOUT) 
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
subroutine H1X_CX(NIN, NOUT, H1, X, XOUT)
  implicit none  
  integer,intent(IN) :: NIN, NOUT                
  double complex, intent(IN)  :: H1(NOUT,NOUT)
  double complex, intent(IN)  :: X(NIN,NOUT)
  double complex, intent(OUT) :: XOUT(NIN,NOUT)

  integer :: i, j

  do i = 1, NOUT
     XOUT(1:NIN, i) = H1(i, 1) * X(1:NIN, 1)
     do j=2, NOUT
         XOUT(1:NIN, i) = XOUT(1:NIN, i) + H1(i, j) * X(1:NIN, j)  
     end do
  end do
end subroutine 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Matrix-vector product used in OSB package:             c
!c                     XOUT = A^T * X                         c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine  H1TX_CX(NIN, NOUT, H1, X, XOUT)
  integer,intent(IN) :: NIN, NOUT                
  double complex,intent(IN) :: H1(NOUT,NOUT)
  double complex,intent(IN) :: X(NIN,NOUT)
  double complex,intent(OUT):: XOUT(NIN,NOUT)

  integer ::    I,J
  
  do i = 1, NOUT
     XOUT(1:NIN, i) = H1(1, i) * X(1:NIN, 1)
     do j=2, NOUT
         XOUT(1:NIN, i) = XOUT(1:NIN, i) + H1(j, i) * X(1:NIN, j)  
     end do
  end do
end subroutine 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Matrix-vector product used in OSB package:             c
!c                     XOUT = A * X                           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine H2X_CX(NIN, NOUT, H1, CORX, X, XOUT)
  implicit none 
  integer,intent(IN) :: NIN, NOUT                
  double complex, intent(IN)  :: H1(NOUT,NOUT)
  double complex, intent(IN)  :: CORX(NIN)
  double complex, intent(IN)  :: X(NIN,NOUT)
  double complex, intent(OUT) :: XOUT(NIN,NOUT)

  integer :: i,j

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
subroutine H2TX_CX(NIN, NOUT, H1, CORX, X, XOUT)
  implicit none 
  integer,intent(IN) :: NIN, NOUT                
  double complex, intent(IN)  :: H1(NOUT,NOUT)
  double complex, intent(IN)  :: CORX(NIN)
  double complex, intent(IN)  :: X(NIN,NOUT)
  double complex, intent(OUT) :: XOUT(NIN,NOUT)

  integer :: i,j

  do i=1,NOUT
     XOUT(1:NIN, i) = H1(1,i) * CORX(1:NIN) * X(1:NIN, 1)
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
subroutine H3X_CX(NIN, NOUT, A, X, XOUT) 
  implicit none
  integer,intent(IN) ::  NIN,NOUT
  double complex,intent(IN) :: A(NIN, NOUT,NOUT)
  double complex,intent(IN)   :: X(NIN, NOUT)
  double complex,intent(OUT)  :: XOUT(NIN, NOUT)

  integer  :: i, j

  do i=1, NOUT     
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
subroutine H3TX_CX(NIN, NOUT, A, X, XOUT) 
  implicit none
  integer,intent(IN) ::  NIN,NOUT
  double complex, intent(IN) :: A(NIN, NOUT, NOUT)
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


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine VTDV_CX(N, D, V, H)
  implicit none
  integer, intent(in)   :: N
  double complex, intent(in)  :: D(N)
  double complex, intent(in)  :: V(N,N)
  double complex, intent(out) :: H(N,N)  

  integer :: i, j

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
subroutine VTHV_CX(N, H0, V, H1)
  implicit none
  integer, intent(in)   :: N
  double complex, intent(in) :: H0(N,N) 
  double complex, intent(in) :: V(N,N)
  double complex, intent(out):: H1(N,N)

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


!cccccccccccccccccccccccccccccccccccccccccccc
!c         Print the Matrix                 c
!cccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine printMat_CX(N, M, mat)
      implicit none
      integer, intent(in) :: N, M
      double complex, intent(in) :: mat(N,M)
    
      integer :: I
      print *, 'Print Complex Matrix: N=', N, ' M=',M

      do I = 1, M
         print *, 'Matrix elements at column=', i
         print *, mat(1:N, i)
      end do

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           B = A^H * A, Matrix-Matrix operation                 c
!c           Used to test GS_FULL_ORTH subroutine                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine AHA_CX(nRow, nCol, A, B)
    implicit none
    integer, intent(IN) :: nRow, nCol
    double complex, intent(IN) :: A(nRow, nCol)
    double complex, intent(OUT):: B(nCol, nCol)  ! B = A^H * A

    integer :: I, J
    
    do I = 1, nCol
       do J = 1, nCol
           B(I, J) = dot_product(A(1:nRow,I), A(1:nRow, J))
       end do
    end do 

end subroutine AHA_CX
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           B = A^T * A, Matrix-Matrix operation              c
!c           Used to test GS_FULL_ORTH subroutine              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine ATA_CX(nRow, nCol, A, B)
    implicit none
    integer, intent(IN) :: nRow, nCol
    double complex, intent(IN) :: A(nRow, nCol)
    double complex, intent(OUT):: B(nCol, nCol)  ! B = A^T * A

    integer :: I, J
    
    do I = 1, nCol
       do J = 1, nCol
           B(I, J) = sum(A(1:nRow,I)* A(1:nRow, J))
       end do
    end do 

end subroutine 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    Test whether a matrix is united orthogonal       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function isOrth_CX(nRow, nCol, A, epsi)
    implicit none
    integer, intent(IN) :: nRow, nCol
    double complex, intent(IN) :: A(nRow, nCol) 
    double precision, intent(IN) :: epsi
    double complex :: tmp

    integer :: I, J

    isOrth_CX = .true.    

    do I = 1, nCol
       do J = 1, nCol         
           tmp = dot_product(A(1:nRow,I), A(1:nRow, J))
           if (I == J) then           ! tmp == 1.0D0 ?
               if (abs(tmp-1.0D0) > EPSI) then
                   isOrth_CX = .false.
                   return
               end if
           else
               if (abs(tmp) > EPSI) then
                   isOrth_CX = .false.
                   return
               end if
           end if
       end do
    end do           

end function 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Maximum Error for I matrix              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function maxMatI_CX(N, A)
    implicit none
    integer, intent(IN) :: N
    double complex, intent(IN) :: A(N, N) 

    integer :: I, J

    maxMatI_CX = 0.0D0
    do i = 1, N
       do J = 1, N
          if (I==J) then
             maxMatI_CX = max(maxMatI_CX, abs(A(I, J) -1.0D0))
          else
             maxMatI_CX = max(maxMatI_CX, abs(A(I, J)))
          end if
       end do
    end do 
end function 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Maximum Error for Symmetric matrix         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function maxSymMat_CX(N, A)
    implicit none
    integer, intent(IN) :: N
    double complex, intent(IN) :: A(N, N)

    integer :: I, J

    maxSymMat_CX = 0.0D0
    do i = 1, N
       do J = I+1, N         
             maxSymMat_CX = max(maxSymMat_CX, abs(A(I, J) -A(J, I)))         
       end do
    end do
  
end function
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Maximum Error for Hermite matrix           c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function maxHerMat_CX(N, A)
    implicit none
    integer, intent(IN) :: N
    double complex, intent(IN) :: A(N, N)

    integer :: I, J

    maxHerMat_CX = 0.0D0
    do i = 1, N
       do J = I+1, N         
             maxHerMat_CX = max(maxHERMat_CX, abs(A(I, J) - conjg(A(J, I))))
       end do
    end do
  
end function
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

