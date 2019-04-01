! Testing program for subroutines in mmt and mtm
program test_mm

!  integer, parameter :: N=500, M=100
  double precision, allocatable, dimension(:,:) :: A, V, U, B, B0, C
  double complex,allocatable,dimension(:,:) :: AX,Vx, Ux, Bx, B0x,Cx
  double precision :: ct1, ct2
  integer :: opt  
  double precision :: isSame, isSameCx
  
  opt=1


  print *, 'Input the size of matrix: N, M'
  read(*,*) N, M

  allocate(A(N,N), V(N,M), U(M,N), B(M,M), B0(M,M), C(N,M), &
           AX(N,N),Vx(N,M), Ux(M,N), Bx(M,M),B0x(M,M),Cx(N,M))

!  select case (opt)
!  case (1)
  print *, ' Testing VTAV for real'  
  call random_number(A); call random_number(V);

  call CPU_Time(ct1)
  call VTAV(N,M,A,V,B)
  call CPU_Time(ct2)
  print *, 'Lapack Time:', (ct2-ct1)

  call CPU_Time(ct1)
  C=matMul(A,V)
  U=transpose(V)
  B0=MatMul(U,C)
  call CPU_Time(ct2)
  print *, 'Intrinsic Time:', (ct2-ct1)

  print *, ' Max. Difference is ',isSame(M,M,B,B0)


  print *, ' Testing VAVT for real'
  call random_number(A); call random_number(U);
  
  call CPU_Time(ct1)
  call VAVT(N,M,A,U,B)
  call CPU_Time(ct2)
  print *, 'Lapack Time:', (ct2-ct1)

  call CPU_Time(ct1)
  V=transpose(U)
  C=matMul(A,V) 
  B0=MatMul(U,C)
  call CPU_Time(ct2)
  print *, 'Intrinsic Time:', (ct2-ct1)

  print *, ' Max Difference is ',isSame(M,M,B,B0)

!  case (2)
  print *, ' Testing VTAV for complex'
  call rand_cx(N,M,Ax);  call rand_cx(N,M,Vx)

  call CPU_Time(ct1)
  call VTAVCx(N,M,Ax,Vx,Bx)
  call CPU_Time(ct2)
  print *, 'Lapack Time:', (ct2-ct1)


  call CPU_Time(ct1)
  Cx=matMul(Ax,Vx)
  Ux=transpose(Vx)
  B0x=MatMul(Ux,Cx)
  call CPU_Time(ct2)
  print *, 'Intrinsic Time:', (ct2-ct1)

  print *, ' Max Difference is ',isSameCx(M,M,Bx,B0x)

  print *, ' Testing VAVT for complex'
  call rand_cx(N,M,Ax);  call rand_cx(M,N,Ux)
  call CPU_Time(ct1)
  call VAVTCx(N,M,Ax,Ux,Bx)
  call CPU_Time(ct2)
  print *, 'Lapack Time:', (ct2-ct1)

  call CPU_Time(ct1)
  Vx=transpose(Ux)
  Cx=matMul(Ax,Vx)
  B0x=MatMul(Ux,Cx)
  call CPU_Time(ct2)
  print *, 'Intrinsic Time:', (ct2-ct1)

  print *, ' Max Difference is ',isSameCx(M,M,Bx,B0x)

!  case (3)
  print *, ' Testing VHAV for complex'
  call rand_cx(N,M,Ax);  call rand_cx(N,M,Vx)

  call CPU_Time(ct1)
  call VHAVCx(N,M,Ax,Vx,Bx)
  call CPU_Time(ct2)
  print *, 'Lapack Time:', (ct2-ct1)

  call CPU_Time(ct1)
  Cx=matMul(Ax,Vx)
  Ux=conjg(transpose(Vx))
  B0x=MatMul(Ux,Cx)
  call CPU_Time(ct2)
  print *, 'Intrinsic Time:', (ct2-ct1)

  print *, ' Max Difference is ',isSameCx(M,M,Bx,B0x)


  print *, ' Testing VTAV for complex'
  call rand_cx(N,M,Ax);  call rand_cx(M,N,Ux)

  call CPU_Time(ct1)
  call VAVHCx(N,M,Ax,Ux,Bx)
  call CPU_Time(ct2)
  print *, 'Lapack Time:', (ct2-ct1)

  call CPU_Time(ct1)
  Vx=Conjg(transpose(Ux))
  Cx=matMul(Ax,Vx)
  B0x=MatMul(Ux,Cx)
  call CPU_Time(ct2)
  print *, 'Intrinsic Time:', (ct2-ct1)

  print *, ' Max Difference is ',isSameCx(M,M,Bx,B0x)


!  end select

  deallocate( A, V, U, B, B0, C, AX, Vx, Ux, Bx, B0x,Cx)

end  


subroutine rand_cx(N, M, A)
  implicit none
  integer, intent(IN) :: N, M
  double complex, intent(OUT) :: A(N,M)

  double precision :: T1(N,M),T2(N,M)

  call random_number(T1)
  call random_number(T2)

  A=CMPLX(T1,T2)
end 

double precision function isSame(N,M,A,B)
   implicit none
   integer, intent(IN) :: N, M
   double precision, intent(IN) :: A(N,M), B(N,M)
   
   integer :: i, j
   double precision :: dab
   isSame=0.0D0

   do i = 1, N
      do j = 1, N
         dab = abs(A(i,j)-B(i,j))
         if (dab>isSame)  isSame=dab        
      end do
   end do

end

double precision function isSameCx(N,M,A,B)
   implicit none
   integer, intent(IN) :: N, M
   double complex, intent(IN) :: A(N,M), B(N,M)
   
   integer :: i, j
   double precision :: dab
   isSameCx=0.0D0

   do i = 1, N
      do j = 1, N
         dab = abs(A(i,j)-B(i,j))
         if (dab>isSameCx) isSameCx=dab         
      end do
   end do

end

