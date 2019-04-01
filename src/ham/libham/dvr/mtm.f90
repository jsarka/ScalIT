!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                B = V^T*A*V                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine VTAV(N, M, A, V, B)   
    implicit none
    integer, intent(IN) :: N, M
    double precision, intent(IN)  :: A(N,N), V(N,M)
    double precision, intent(OUT) :: B(M,M)

    double precision, allocatable :: tmp(:,:)

    allocate(tmp(N, M))
 
    call DGEMM('N','N',N,M,N,1.0D0,A,N,V,N,0.0D0,tmp,N)

    call DGEMM('T','N',M,M,N,1.0D0,V,N,tmp,N,0.0D0,B,M)

    deallocate(tmp)

end
!************************************
subroutine VTAV0(N, A, V, B)   
    implicit none
    integer, intent(IN) :: N
    double precision, intent(IN)  :: A(N,N), V(N,N)
    double precision, intent(OUT) :: B(N,N)

    call VTAV(N, N, A, V, B)
end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!************************************
subroutine VTAVCx(N, M, ACx, VCx, BCx)   
    implicit none
    integer, intent(IN) :: N, M
    double complex, intent(IN)  :: ACx(N,N), VCx(N,M)
    double complex, intent(OUT) :: BCx(M,M)

    double complex, allocatable :: tmp(:,:)
    double complex :: alpha=(1.0D0,0.0D0), beta=(0.0D0,0.0D0)

    allocate(tmp(N, M))
 
    call ZGEMM('N','N',N,M,N,alpha,ACx,N,VCx,N,beta,tmp,N)

    call ZGEMM('T','N',M,M,N,alpha,VCx,N,tmp,N,beta,BCx,M)

    deallocate(tmp)
end
!************************************
subroutine VTAV0Cx(N, ACx, VCx, BCx)   
    implicit none
    integer, intent(IN) :: N
    double complex, intent(IN)  :: ACx(N,N), VCx(N,N)
    double complex, intent(OUT) :: BCx(N,N)

    call VTAVCx(N, N, ACx, VCx, BCx)
end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!************************************
subroutine VHAVCx(N, M, ACx, VCx, BCx)   
    implicit none
    integer, intent(IN) :: N, M
    double complex, intent(IN)  :: ACx(N,N), VCx(N,M)
    double complex, intent(OUT) :: BCx(M,M)

    double complex, allocatable :: tmp(:,:)
    double complex :: alpha=(1.0D0,0.0D0), beta=(0.0D0,0.0D0)

    allocate(tmp(N, M))  

    call ZGEMM('N','N',N,M,N,alpha,ACx,N,VCx,N,beta,tmp,N)

    call ZGEMM('C','N',M,M,N,alpha,VCx,N,tmp,N,beta,BCx,M)

    deallocate(tmp)
end
!************************************
subroutine VHAV0Cx(N, ACx, VCx, BCx)   
    implicit none
    integer, intent(IN) :: N
    double complex, intent(IN)  :: ACx(N,N), VCx(N,N)
    double complex, intent(OUT) :: BCx(N,N)

    call VHAVCx(N, N, ACx, VCx, BCx)
end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     H = V^T * D * V,  D is a diagonal matrix      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine VTDV(N, M, D, V, B)
  implicit none
  integer, intent(IN) :: N,M
  double precision, intent(IN) :: D(N), V(N,M)
  double precision, intent(OUT):: B(M,M) 

  integer :: i, j

  do i = 1, M
     do j = 1, I
        B(I, J) = Sum(V(1:N,I)*D(1:N)*V(1:N, J))
        B(J, I) = B(I, J)
     end do
  end do        

end
!************************************
subroutine VTDV0(N, D, V, B)
  implicit none
  integer, intent(IN) :: N
  double precision, intent(IN) :: D(N), V(N,N)
  double precision, intent(OUT):: B(N,N)

  call VTDV(N, N, D, V, B) 

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine VTDVCx(N, M, D, VCx, BCx)
  implicit none
  integer, intent(IN) :: N,M
  double precision, intent(IN) :: D(N)
  double complex, intent(IN)   :: VCx(N,M)
  double complex, intent(OUT)  :: BCx(M,M) 

  integer :: i, j

  do i = 1, M
     do j = 1, M
        BCx(I, J) = Sum(VCx(1:N,I)*D(1:N)*VCx(1:N, J))    
     end do
  end do        

end
!************************************
subroutine VTDV0Cx(N, D, VCx, BCx)
  implicit none
  integer, intent(IN) :: N
  double precision, intent(IN) :: D(N)
  double complex, intent(IN)   :: VCx(N,N)
  double complex, intent(OUT)  :: BCx(N,N)

  call VTDVCx(N, N, D, VCx, BCx) 

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine VHDVCx(N, M, D, VCx, BCx)
  implicit none
  integer, intent(IN) :: N,M
  double precision, intent(IN) :: D(N)
  double complex, intent(IN)   :: VCx(N,M)
  double complex, intent(OUT)  :: BCx(M,M) 

  integer :: i, j

  do i = 1, M
     do j = 1, I
        BCx(I, J) = Sum(Conjg(VCx(1:N,I))*D(1:N)*VCx(1:N, J))
        BCx(J, I) = Conjg(BCx(I, J))
     end do
  end do        

end
!************************************
subroutine VHDV0Cx(N, D, VCx, BCx)
  implicit none
  integer, intent(IN) :: N
  double precision, intent(IN) :: D(N)
  double complex, intent(IN)   :: VCx(N,N)
  double complex, intent(OUT)  :: BCx(N,N)

  call VHDVCx(N, N, D, VCx, BCx) 

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc



