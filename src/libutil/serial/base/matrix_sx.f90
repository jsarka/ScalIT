!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine VTDV_HX(N, D, V, H)
  implicit none
  integer, intent(in)   :: N
  double complex, intent(in)  :: D(N)
  double complex, intent(in)  :: V(N,N)
  double complex, intent(out) :: H(N,N)  

  integer :: i, j

  do i = 1, N
     do j = 1, N
        H(i, j) = sum( CONJG((V(1:N,i))) * D(1:N) * V(1:N,j) )
     end do
  end do        
end subroutine 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Calculate H1 = V^H * H0 * V,  D is a normal matrix    c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   H(i, j) = sum(k) [V(k,i) * H0(k, l) * V(k, j)]       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine VTHV_HX(N, H0, V, H1)
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
           tmp = tmp + sum( CONJG(V(1:N, i)) * H0(1:N, k)) * V(k, j)
        end do
        H1(i, j) = tmp
     end do
  end do        
end subroutine 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    Test whether a matrix is united orthogonal       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function isOrth_SX (nRow, nCol, A, epsi)
    implicit none
    integer, intent(IN) :: nRow, nCol
    double complex, intent(IN) :: A(nRow, nCol) 
    double precision, intent(IN) :: epsi
    double complex :: tmp

    integer :: I, J
    double complex :: dot_cx

    isOrth_SX = .true.    

    do I = 1, nCol
       do J = 1, nCol         
           tmp = sum(A(1:nRow,I)*A(1:nRow, J))
           if (I == J) then           ! tmp == 1.0D0 ?
               if (abs(tmp-1.0D0) > EPSI) then
                   isOrth_SX = .false.
                   return
               end if
           else
               if (abs(tmp) > EPSI) then
                   isOrth_SX = .false.
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
double precision function maxMatI_SX(N, A)
    implicit none
    integer, intent(IN) :: N
    double complex, intent(IN) :: A(N, N) 

    integer :: I, J

    maxMatI_SX = 0.0D0
    do i = 1, N
       do J = 1, N
          if (I==J) then
             maxMatI_SX = max(maxMatI_SX, abs(A(I, J) -1.0D0))
          else
             maxMatI_SX = max(maxMatI_SX, abs(A(I, J)))
          end if
       end do
    end do 
end function 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

