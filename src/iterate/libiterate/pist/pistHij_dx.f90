!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Create H Matrix in PIST                  c
!c           H(I, J) = <VI|H|VJ>                      c
!c              Complex Version                       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function pistHij_dx(N, X, M, H0XDX, LinSolv, HMAT)
   implicit none
   integer, intent(IN)   :: N, M
   double precision, intent(IN)    :: X(N)
   external             :: H0XDX
   integer, external    :: LinSolv
   double complex, intent(OUT) :: HMAT(M,M)

   double precision :: VJ(N, M), WJ(N)
   double precision :: TMP, RES
   double complex   :: WJ0(N), WJ1(N)
   integer :: I, J, NUM
   logical :: MGS_ORTH
  
   tmp  = DSQRT(dot_product(X(1:N), X(1:N)))

   pistHIJ_DX = 0

   if (tmp == 0.0D0)    return

   VJ(1:N, 1)   = X(1:N)/tmp  

   PISTHIJ_DX = M

   do J = 1, M-1          ! WJ = VJ(1:N, J) / (H-EI)

      NUM = LinSolv(N, VJ(1:N, J), WJ, RES)          
                        
      if (NUM <= 0) then       
          PISTHIJ_DX = num;      return
      end if       

      if (.not. MGS_ORTH(N, WJ, J, VJ)) then        
          PISTHIJ_DX = -J;       return        
      end if    ! VJ(1:N, J+1) = H*VJ

   end do
 
   do J = 1, M
      call H0XDX(N, VJ(1:N, J), WJ0)       
      do I = 1, M
         WJ1(1:N) = VJ(1:N,I)*WJ0(1:N)
         HMAT(I, J)  = sum(WJ1(1:N))  
      end do
   end do

end function PISTHIJ_DX
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Create H Matrix in PIST                  c
!c           H(I, J) = <VI|H|VJ>                      c
!c            Complex Version                         c          
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function pistH0_DX(N, X, M_MAX, M_LOW, M_HIGH, H0XDX, LinSolv, VJ, HMAT)
   implicit none
   integer, intent(IN)   :: N,  M_MAX, M_LOW, M_HIGH
   double precision, intent(IN)    :: X(N)
   external             :: H0XDX
   integer, external    :: LinSolv
   double precision, intent(INOUT) :: VJ(N, M_HIGH)
   double complex, intent(INOUT) :: HMAT(M_MAX, M_HIGH)

   double precision :: WJ(N)
   double complex   :: WJ0(N), WJ1(N)
   double precision :: TMP, RES
   integer :: I, J, NUM, M_START
   logical :: MGS_ORTH
  
   pistH0_DX= 0

   if ((M_LOW<1) .OR. (M_LOW>=M_HIGH) .OR. (M_HIGH>M_MAX) )     return

   if (M_LOW == 1) then
      M_START = 1
      tmp  = DSQRT(dot_product(X(1:N), X(1:N)))

      if (tmp == 0.0D0)  return

      VJ(1:N, 1)   = X(1:N)/tmp  
   else
      M_START = M_LOW - 1
   end if

   PISTH0_DX = M_HIGH
   do J = M_START, M_HIGH-1          ! WJ = VJ(1:N, J) / (H-EI)
      NUM = LinSolv( N, VJ(1:N, J), WJ, RES)
                                       
      if (NUM <= 0) then       
         PISTH0_DX = num;     return
      end if       

      if (.not. MGS_ORTH(N, WJ, J, VJ)) then        
          PISTH0_DX = -J;      return        
      end if    ! VJ(1:N, J+1) = H*VJ
   end do
 
   do J = M_LOW, M_HIGH
      call H0XDX(N, VJ(1:N, J), WJ0)    
      do I = 1, M_HIGH
         WJ1(1:N) = VJ(1:N, I)*WJ0(1:N)
         HMAT(I, J)  = sum(WJ1(1:N))  
      end do
   end do

   do J = 1, M_LOW-1 
      call H0XDX(N, VJ(1:N, J), WJ0)  
      do I = M_LOW, M_HIGH
         WJ1(1:N) = VJ(1:N, I)*WJ0(1:N)
         HMAT(I, J)  = sum(WJ1(1:N))  
      end do
   end do

end function PISTH0_DX
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


