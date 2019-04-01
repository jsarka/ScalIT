!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                    makeI                               c
!c   Create a Identity matrix: A(i,i)=1, A(i,j) = 0       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine makeI(N, A)
  implicit none
  integer, intent(IN) :: N
  double precision, intent(OUT) :: A(N, N)

  integer :: i

  A(1:N, 1:N) = 0.0D0
  do  i = 1, N
     A(i, i)   = 1.0D0
  end do
end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccc
!c         Randomly generated Matrix           c
!ccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine randMat(n,m, mat)
    implicit none
    integer, intent(IN) :: N, M
    double precision, intent(OUT) :: mat(N,M)    

    call random_number(mat)
end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         Fill the Diagonal elements of matrix      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine fillDiagMat(NIN, NOUT, diagMat, Mat)
     implicit none
     integer, intent(IN) :: NIN, NOUT
     double precision, intent(IN) :: diagMat(NIN, NOUT, NOUT) 
     double precision, intent(OUT):: mat(NIN*NOUT, NIN*NOUT) 

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
subroutine addDiagMat(NIN, NOUT, diagMat, Mat)
     implicit none
     integer, intent(IN) :: NIN, NOUT
     double precision, intent(IN)    :: diagMat(NIN, NOUT, NOUT)
     double precision, intent(INOUT) :: mat(NIN*NOUT, NIN*NOUT) 

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
subroutine H1X(NIN, NOUT, H1, X, XOUT)
  implicit none  
  integer,intent(IN) :: NIN, NOUT                
  double precision, intent(IN) :: H1(NOUT,NOUT)
  double precision, intent(IN) :: X(NIN,NOUT)
  double precision, intent(OUT) :: Xout(NIN,NOUT)

  integer :: i,j

  do i = 1, NOUT
     XOUT(1:NIN, i) = H1(i, 1) * X(1:NIN, 1)
     do j=2,NOUT
         XOUT(1:NIN, i) = XOUT(1:NIN, i) + H1(i, j) * X(1:NIN, j)  
     end do
  end do

end 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Matrix-vector product used in OSB package:             c
!c                     XOUT = A^T * X                         c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine  H1TX(NIN, NOUT, H1, X, XOUT)
  integer,intent(IN) :: NIN, NOUT                
  double precision, intent(IN) :: H1(NOUT,NOUT)
  double precision, intent(IN) :: X(NIN,NOUT)
  double precision, intent(OUT) :: Xout(NIN,NOUT)

  integer ::    I, J

  do i = 1, NOUT
     XOUT(1:NIN, i) = H1(1, i) * X(1:NIN, 1)
     do j=2,NOUT
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
subroutine H2X(NIN, NOUT, H1, CORX, X, XOUT)
  implicit none 
  integer,intent(IN) :: NIN, NOUT                
  double precision, intent(IN)  :: H1(NOUT,NOUT)
  double precision, intent(IN)  :: CORX(NIN)
  double precision, intent(IN)  :: X(NIN,NOUT)
  double precision, intent(OUT) :: Xout(NIN,NOUT)

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
subroutine H2TX(NIN, NOUT, H1, CORX, X, XOUT)
  implicit none 
  integer,intent(IN) :: NIN, NOUT                
  double precision, intent(IN)  :: H1(NOUT,NOUT)
  double precision, intent(IN)  :: CORX(NIN)
  double precision, intent(IN)  :: X(NIN,NOUT)
  double precision, intent(OUT) :: Xout(NIN,NOUT)

  integer :: i,j

  do i=1,NOUT
     XOUT(1:NIN, i) = H1(1,i) * CORX(1:NIN) * X(1:NIN, 1)
     do j=2,NOUT
         XOUT(1:NIN,i) = XOUT(1:NIN,i)+H1(j,i)*CORX(1:NIN)*X(1:NIN,j)
     end do
  end do
end subroutine 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     matrix-vector product:  XOUT = A^T * X       c
!c     XOUT, X must be different                    c
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine H3X(NIN, NOUT, A, X, XOUT) 
  integer,intent(IN) ::  NIN,NOUT
  double precision, intent(IN)  :: A(NIN, NOUT, NOUT)
  double precision, intent(IN)  :: X(NIN, NOUT)
  double precision, intent(OUT) :: Xout(NIN, NOUT)

  integer  :: i, j
      
  do i=1, NOUT        
     XOUT(1:NIN, i) = A(1:NIN, i, 1) * X(1:NIN, 1)
     do j=2, NOUT
          XOUT(1:NIN, i) = XOUT(1:NIN, i) + A(1:NIN, i, j) * X(1:NIN, j)
     end do
  end do  

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     matrix-vector product:  XOUT = A^T * X       c
!c     XOUT, X must be different                    c
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine H3TX(NIN, NOUT, A, X, XOUT) 
  integer,intent(IN) ::  NIN,NOUT
  double precision, intent(IN)  :: A(NIN, NOUT, NOUT)
  double precision, intent(IN)  :: X(NIN, NOUT)
  double precision, intent(OUT) :: Xout(NIN, NOUT)

  integer  :: j, k      

  do i=1, NOUT 
     XOUT(1:NIN, i) = A(1:NIN, 1, i) * X(1:NIN, 1)
     do j=2, NOUT
          XOUT(1:NIN, i) = XOUT(1:NIN, i) + A(1:NIN, j, i) * X(1:NIN, j)
     end do
  end do  
end subroutine 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccc
!c         Print the Matrix                 c
!cccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine printMat(N, M, mat)
      implicit none
      integer, intent(in) :: N, M
      double precision, intent(in) :: Mat(N, M)
    
      integer :: I
      print *, 'Print Real Matrix: N=', N, ' M=',M

      do I = 1, M
         print *, 'Matrix elements at column=', i
         print *, mat(1:N, i)
      end do

      end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           B = A^T * A, Matrix-Matrix operation               c
!c           Used to test GS_FULL_ORTH subroutine               c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine ATA (nRow, nCol, A, B)
    implicit none
    integer, intent(IN) :: nRow, nCol
    double precision, intent(IN)  :: A(nRow, nCol)
    double precision, intent(OUT) :: B(nCol, nCol)

    integer :: I, J
    
    do I = 1, nCol
       do J = 1, nCol
           B(I, J) = dot_product(A(1:nRow,I), A(1:nRow, J))
       end do
    end do 
end subroutine ATA
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    Test whether a matrix is united orthogonal       v
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function isOrth (nRow, nCol, A, epsi)
    implicit none
    integer, intent(IN) :: nRow, nCol
    double precision, intent(IN) :: A(nRow, nCol)
    double precision, intent(IN) :: epsi
    double precision :: tmp

    integer :: I, J

    isOrth = .true.    

    do I = 1, nCol
       do J = 1, nCol         
           tmp = dot_product(A(1:nRow,I), A(1:nRow, J))
           if (I == J) then           ! tmp == 1.0D0 ?
               if (abs(tmp-1.0D0) > EPSI) then
                   isOrth = .false.
                   return
               end if
           else
               if (abs(tmp) > EPSI) then
                   isOrth = .false.
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
double precision function maxMatI(N, A)
    implicit none
    integer, intent(IN) :: N
    double precision, intent(IN) :: A(N, N)

    integer :: I, J

    maxMatI = 0.0D0
    do i = 1, N
       do J = 1, N
          if (I==J) then
             maxMatI = max(maxMatI, abs(A(I, J) -1.0D0))
          else
             maxMatI = max(maxMatI, abs(A(I, J)))
          end if
       end do
    end do  
end function
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Maximum Error for Symmetric matrix         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function maxSymMat(N, A)
    implicit none
    integer, intent(IN) :: N
    double precision, intent(IN) :: A(N, N)

    integer :: I, J

    maxSymMat = 0.0D0
    do i = 1, N
       do J = I+1, N         
             maxSymMat = max(maxSymMat, abs(A(I, J) -A(J, I)))         
       end do
    end do  
end function
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




