!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      Create H Matrix in PIST: H(I,J) = <VI|H|VJ>        c
!c      VI is the PIST vectors:                            c
!c           WJ = (H-E)^-1*VI(I)                           c
!c           VI(I+1) = WJ - Sum((VI(K):WJ)*VI(K))          c
!c                  k = 1, 2, ..., I                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function pistHij(N, X, M, H0X, LinSolv, HMAT)
   implicit none
   integer, intent(IN)   :: N, M
   double precision, intent(IN)    :: X(N)
   external           :: H0X
   integer, external  :: LinSolv     ! Linear Solver
   double precision, intent(OUT) :: HMAT(M,M)

   double precision :: VJ(N,M),WJ(N)
   double precision :: TMP, RES
   integer :: I, J, NUM
   logical :: MGS_ORTH
  
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   tmp  = DSQRT(dot_product(X(1:N), X(1:N)))

   pistHIJ = 0

   if (tmp == 0.0D0)        return

   VJ(1:N, 1)   = X(1:N)/tmp  

   PISTHIJ = M

   do J = 1, M-1          ! WJ = VJ(1:N, J) / (H-EI)

      NUM = LinSolv( N, VJ(1:N, J), WJ, RES)   

      if (NUM <= 0) then       
         PISTHIJ = num;  return
      end if       

      if (.not. MGS_ORTH(N, WJ, J, VJ)) then        
         PISTHIJ = -J;      return        
      end if      ! VJ(1:N, J+1) = H*VJ

   end do
 
   do J = 1, M
      call H0X(N, VJ(1:N, J), WJ)  ! WJ = H * VJ
      do I = J, M                  ! H(I,J)=<UI|H|UJ>
        HMAT(I, J)  = dot_product(VJ(1:N,I), WJ)  
        HMAT(J, I)  = HMAT(I, J)
      end do
   end do

end function PISTHIJ
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Create H Matrix in PIST: H(I, J) = <VI|H|VJ>        c
!c       VI is the PIST vectors:                           c   
!c           WJ = (H-E)^-1*VI(I)                           c
!c           VI(I+1) = WJ - Sum((VI(K):WJ)*VI(K))          c
!c                    k = 1, 2, ..., I                     c
!c   This version need to store the old vectors and HMAT,  c
!c   but it will save the time to calculate the previous   c
!c   calculated HIJ and vectors. The previous calculated   c
!c   ones are HMAT(1:M_LOW-1, 1:M_LOW-1), VJ(N, 1:M_LOW-1) c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function pistH0(N, X, M_MAX, M_LOW, M_HIGH, H0X, LinSolv, VJ, HMAT)
   implicit none
   integer, intent(IN)   :: N, M_MAX, M_LOW, M_HIGH
   double precision, intent(IN)    :: X(N)
   external          :: H0X
   integer, external :: LinSolv
   double precision, intent(INOUT) :: VJ(N, M_HIGH), HMat(M_MAX, M_HIGH)

   double precision :: WJ(N)
   double precision :: TMP, RES
   integer :: I, J, NUM, M_START
   logical :: MGS_ORTH
  
   pistH0= 0

   if ((M_LOW<1) .OR. (M_LOW>=M_HIGH) .OR. (M_HIGH>M_MAX) )   return

   if (M_LOW == 1) then
      M_START = 1

      tmp  = DSQRT(dot_product(X(1:N), X(1:N)))

      if (tmp == 0.0D0)  return

      VJ(1:N, 1)   = X(1:N)/tmp  
   else
      M_START = M_LOW - 1
   end if

   PISTH0 = M_HIGH
  
   do J = M_START, M_HIGH-1          ! WJ = VJ(1:N, J) / (H-EI)

      num = LinSolv(N, VJ(1:N, J), WJ, RES)

      print *,'# of Lanczos iter.:',J, '. # of iter. for Linear Solver:', num
    
      if (NUM <= 0) then       
         PISTH0 = num;   return
      end if       

      if (.not. MGS_ORTH(N, WJ, J, VJ)) then        
          PISTH0 = -J;       return        
      end if    ! VJ(1:N, J+1) = H*VJ

   end do
 
   do J = M_LOW, M_HIGH
      call H0X(N, VJ(1:N, J), WJ)     ! WJ = H * VJ
      do I = 1, M_HIGH                ! H(I,J)=<UI|H|UJ>
        HMAT(I, J)  = dot_product(VJ(1:N,I), WJ)       
        HMAT(J, I)  = HMAT(I, J)
      end do
   end do

end function PISTH0
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



