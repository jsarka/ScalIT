!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                B = V*A*V^T                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine VAVT(N, M, A, V, B)   
    implicit none
    integer, intent(IN) :: N, M
    double precision, intent(IN)  :: A(N,N), V(M,N)
    double precision, intent(OUT) :: B(M,M)

    double precision, allocatable :: tmp(:,:)

    allocate(tmp(N, M))
 
    call DGEMM('N','T',N,M,N,1.0D0,A,N,V,M,0.0D0,tmp,N)

    call DGEMM('N','N',M,M,N,1.0D0,V,M,tmp,N,0.0D0,B,M)

    deallocate(tmp)

end
!************************************
subroutine VAVT0(N, A, V, B)   
    implicit none
    integer, intent(IN) :: N
    double precision, intent(IN)  :: A(N,N), V(N,N)
    double precision, intent(OUT) :: B(N,N)

    call VAVT(N, N, A, V, B)
end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!************************************
subroutine VAVTCx(N, M, ACx, VCx, BCx)   
    implicit none
    integer, intent(IN) :: N, M
    double complex, intent(IN)  :: ACx(N,N), VCx(M,N)
    double complex, intent(OUT) :: BCx(M,M)

    double complex :: OneCx=(1.0D0,0.0D0),ZeroCx=(0.0D0,0.0D0)
    double complex, allocatable :: tmp(:,:)

    allocate(tmp(N, M))
 
    call ZGEMM('N','T',N,M,N,OneCx,ACx,N,VCx,M,ZeroCx,tmp,N)

    call ZGEMM('N','N',M,M,N,OneCx,VCx,M,tmp,N,ZeroCx,BCx,M)

    deallocate(tmp)
end
!************************************
subroutine VAVT0Cx(N, ACx, VCx, BCx)   
    implicit none
    integer, intent(IN) :: N
    double complex, intent(IN)  :: ACx(N,N), VCx(N,N)
    double complex, intent(OUT) :: BCx(N,N)

    call VAVTCx(N, N, ACx, VCx, BCx)
end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!************************************
subroutine VAVHCx(N, M, ACx, VCx, BCx)   
    implicit none
    integer, intent(IN) :: N, M
    double complex, intent(IN)  :: ACx(N,N), VCx(M,N)
    double complex, intent(OUT) :: BCx(M,M)

    double complex,parameter :: OneCx=(1.0D0,0.0D0), ZeroCx=(0.0D0,0.0D0)
    double complex, allocatable :: tmp(:,:)

    allocate(tmp(N, M))
 
    call ZGEMM('N','C',N,M,N,OneCx,ACx,N,VCx,M,ZeroCx,tmp,N)

    call ZGEMM('N','N',M,M,N,OneCx,VCx,M,tmp,N,ZeroCx,BCx,M)

    deallocate(tmp)
end
!************************************
subroutine VAVH0Cx(N, ACx, VCx, BCx)   
    implicit none
    integer, intent(IN) :: N
    double complex, intent(IN)  :: ACx(N,N), VCx(N,N)
    double complex, intent(OUT) :: BCx(N,N)

    call VAVHCx(N, N, ACx, VCx, BCx)
end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     H = V^T * D * V,  D is a diagonal matrix      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine VDVT(N, M, D, V, B)
  implicit none
  integer, intent(IN) :: N,M
  double precision, intent(IN) :: D(N), V(M,N)
  double precision, intent(OUT):: B(M,M) 

  integer :: i, j

  do i = 1, M
     do j = 1, I
        B(I, J) = Sum(V(I,1:N)*D(1:N)*V(J,1:N))
        B(J, I) = B(I, J)
     end do
  end do        

end
!************************************
subroutine VDVT0(N, D, V, B)
  implicit none
  integer, intent(IN) :: N
  double precision, intent(IN) :: D(N), V(N,N)
  double precision, intent(OUT):: B(N,N)

  call VDVT(N, N, D, V, B) 

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine VDVTCx(N, M, D, VCx, BCx)
  implicit none
  integer, intent(IN) :: N,M
  double precision, intent(IN) :: D(N)
  double complex, intent(IN)   :: VCx(M,N)
  double complex, intent(OUT)  :: BCx(M,M) 

  integer :: i, j

  do i = 1, M
     do j = 1, M
        BCx(I, J) = Sum(VCx(I,1:N)*D(1:N)*VCx(J,1:N))    
     end do
  end do        

end
!************************************
subroutine VDVT0Cx(N, D, VCx, BCx)
  implicit none
  integer, intent(IN) :: N
  double precision, intent(IN) :: D(N)
  double complex, intent(IN)   :: VCx(N,N)
  double complex, intent(OUT)  :: BCx(N,N)

  call VDVTCx(N, N, D, VCx, BCx) 

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine VDVHCx(N, M, D, VCx, BCx)
  implicit none
  integer, intent(IN) :: N,M
  double precision, intent(IN) :: D(N)
  double complex, intent(IN)   :: VCx(M,N)
  double complex, intent(OUT)  :: BCx(M,M) 

  integer :: i, j

  do i = 1, M
     do j = 1, I
        BCx(I, J) = Sum(VCx(I,1:N)*D(1:N)*Conjg(VCx(J,1:N)))
        BCx(J, I) = Conjg(BCx(I, J))
     end do
  end do        

end
!************************************
subroutine VDVH0Cx(N, D, VCx, BCx)
  implicit none
  integer, intent(IN) :: N
  double precision, intent(IN) :: D(N)
  double complex, intent(IN)   :: VCx(N,N)
  double complex, intent(OUT)  :: BCx(N,N)

  call VDVHCx(N, N, D, VCx, BCx) 

end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc



