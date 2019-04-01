!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    Subroutines to calculate the transform matrix        c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function getDVRT1(fname1, fname2, N, M, T)
   character(len=*), intent(IN) :: fname1, fname2
   integer, intent(IN)  :: N, M
   double precision, intent(OUT) :: T(N, M)

   integer :: info, num
   double precision, allocatable :: U(:,:), V(:,:), W(:,:)
   double precision, allocatable :: X10(:),Heig(:),X20(:),H2(:,:) 
   logical :: dvr_step30

   getDVRT1 = .false. 
   if ( N < M ) return 

   open(99, file=fname2, status='old', form='unformatted', iostat=info)
   if (info/=0) return
   read(99) num
   if (num /= N) then
      close(99); return
   end if

   allocate (V(N, M), W(M,M),X10(N),HEig(N),X20(M),H2(M,M), stat=info)
   if (info/=0) then
      close(99); return
   end if

   read(99) X10(1:N), Heig(1:N), V(1:N,1:M)
   close(99)

   if (.NOT.DVR_Step30(N, M, X10, V, HEig, X20, H2, W)) then
       deallocate(V, W, X10, HEig,X20, H2);   return
   end if
    
   deallocate(HEig,X20,H2)
   T(1:N,1:M)=0.0D0      ! T(N,M) = V(N,M)*W(M,M)
   call dgemm('N','N',N,M,M,1.0D0,V,N,W,M,0.0D0,T,N)
   V(1:N,1:M) = T(1:N,1:M)

  deallocate(W)
  open(99, file=fname1, status='old', form='unformatted', iostat=info)
  if (info/=0) then
      deallocate(V, X10); return
  end if
  read(99) num
  if (num /= N) then
     close(99); deallocate(V,X10); return
  end if

  allocate(U(N,N),H2(N,N), stat=info)
  if (info/=0) then
     close(99); deallocate(V,X10); return
  end if

  read(99) X10(1:N), H2(1:N,1:N), U(1:N,1:N)
 
  close(99)

  deallocate(X10, H2)    ! T(N,M) = U(N,N)*V(N,M) 
  call dgemm('N','N',N,M,N,1.0D0,U,N,V,N,0.0D0,T,N)

  deallocate(U, V)

  getDVRT1 = .false.

end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function getDVRT2(fname, N, M, T)
   character(len=*),intent(IN) :: fname
   integer, intent(IN)  :: N, M 
   double precision, intent(OUT) :: T(N, M)
   
   integer :: info, num
   double precision, allocatable :: V(:,:), W(:,:)
   double precision, allocatable :: X10(:),Heig(:),X20(:),H2(:,:)
   logical :: dvr_step30   

   getDVRT2 = .false.
   if ( N < M ) return
   
   open(99, file=fname, status='old', form='unformatted', iostat=info)
   if (info/=0) return
   read(99) num 
   if (num /= N) then
      close(99); return
   end if
   
   allocate (V(N, M), W(M,M),X10(N),HEig(N),X20(M),H2(M,M), stat=info)
   if (info/=0) then
      close(99); return
   end if
   
   read(99) X10(1:N), Heig(1:N), V(1:N,1:M)
   close(99)
   
   if (.NOT.DVR_Step30(N, M, X10, V, HEig, X20, H2, W)) then
       deallocate(V, W, X10, HEig,X20, H2);   return
   end if
    
   T(1:N,1:M)=0.0D0    ! T(N,M)=V(N,M)*W(M,M)
   call dgemm('N','N',N,M,N,1.0D0,V,N,W,M,0.0D0,T,N)
   
   deallocate(HEig, X20, H2, V) 

end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


