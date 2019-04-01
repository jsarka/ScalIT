!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Create H Matrix in PIST                  c
!c           H(I, J) = <VI*|H|VJ>                      c
!c              Complex Version                       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function pistHij_sx(N, X, M, H0X, LinSolv, HMAT)
   implicit none
   integer, intent(IN)   :: N, M
   double complex, intent(IN) :: X(N)
   external             :: H0X
   integer, external    :: LinSolv
   double complex, intent(OUT) :: HMAT(M,M)

   double complex   :: VJ(N, M), WJ(N), dot_cx
   double complex   :: tmp
   double precision :: RES
   integer :: I, J, NUM
   logical :: MGS_ORTH_SX
  
   tmp  = SQRT(dot_cx(N,X(1:N), X(1:N)))

   pistHIJ_SX = 0

   if (tmp == 0.0D0)    return

   VJ(1:N, 1)   = X(1:N)/tmp  

   PISTHIJ_SX = M

   do J = 1, M-1          ! WJ = VJ(1:N, J) / (H-EI)

      NUM = LinSolv(N, VJ(1:N, J), WJ, RES)          
                        
      if (NUM <= 0) then       
         PISTHIJ_SX = num;     return
      end if       

      if (.not. MGS_ORTH_SX(N, WJ, J, VJ)) then        
         PISTHIJ_SX = -J;       return        
      end if    ! VJ(1:N, J+1) = H*VJ

   end do
 
   do J = 1, M
      call H0X(N, VJ(1:N, J), WJ)         ! WJ = H * VJ
      do I = 1, M
        HMAT(I, J)  = dot_cx(N, VJ(1:N,I), WJ(1:N)) 
      end do
   end do

end function PISTHIJ_SX
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Create H Matrix in PIST                  c
!c           H(I, J) = <VI|H|VJ>                      c
!c            Complex Version                         c          
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function pistH0_SX(N, X, M_MAX, M_LOW, M_HIGH, H0X, LinSolv, VJ, HMAT)
   implicit none
   integer, intent(IN)   :: N,  M_MAX, M_LOW, M_HIGH
   double complex, intent(IN)    :: X(N)
   external             :: H0X
   integer, external    :: LinSolv
   double complex, intent(INOUT) :: VJ(N, M_HIGH), HMAT(M_MAX, M_HIGH)

   double complex   :: WJ(N)
   double complex   :: dot_cx, tmp
   double precision :: RES
   integer :: I, J, NUM, M_START
   logical :: MGS_ORTH_SX
  
   pistH0_SX= 0

   if ((M_LOW<1) .or. (M_LOW>=M_HIGH) .or. (M_HIGH>M_MAX) )     return

   if (M_LOW == 1) then
      M_START = 1
      tmp  = SQRT(dot_cx(N,X(1:N), X(1:N)))
     
      if (tmp == 0.0D0)    return
     
      VJ(1:N, 1)   = X(1:N)/tmp  
   else
      M_START = M_LOW - 1
   end if

   PISTH0_SX = M_HIGH
   do J = M_START, M_HIGH-1          ! WJ = VJ(1:N, J) / (H-EI)
      NUM = LinSolv( N, VJ(1:N, J), WJ, RES)
                                       
      if (NUM <= 0) then       
         PISTH0_SX = num;     return
      end if       

      if (.not. MGS_ORTH_SX(N, WJ, J, VJ)) then        
         PISTH0_SX = -J;       return        
      end if    ! VJ(1:N, J+1) = H*VJ
   end do
 
   do J = M_LOW, M_HIGH
      call H0X(N, VJ(1:N, J), WJ)         ! WJ = H * VJ
      do I = 1, M_HIGH
         HMAT(I, J)  = dot_cx(N,VJ(1:N,I), WJ)  ! HMAT(I,J) = <UI|H|UJ>        
      end do
   end do

   do J = 1, M_LOW-1 
      call H0X(N, VJ(1:N, J), WJ)  ! VJ = H * U
      do I = M_LOW, M_HIGH
          HMAT(I, J)  = dot_cx(N,VJ(1:N,I), WJ(1:N))  ! HMAT(I,J) = <UI|H|UJ>
      end do
   end do

end function PISTH0_SX
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


