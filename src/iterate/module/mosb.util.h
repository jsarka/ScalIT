!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             Some Utility subroutines                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!*************************************************************
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine randVec_cx(N, X)
   integer, intent(IN) :: N
   double complex, intent(OUT) :: X(N)

   double precision :: A(N), B(N)

   call random_number(A)

   call random_number(B)

   X(1:N) = DCMPLX(A(1:N), B(1:N))

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!*************************************************************
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     subroutine copyVec(N, X0, X1)
        integer, intent(IN) :: N
        double precision, intent(IN) :: X0(N)
        double precision, intent(OUT):: X1(N)
        
        X1(1:N) = X0(1:N)

     end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     subroutine copyVecCx(N, X0, X1)
        integer, intent(IN) :: N
        double complex, intent(IN) :: X0(N)
        double complex, intent(OUT):: X1(N)
        
        X1(1:N) = X0(1:N)

     end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine copy2LVec(N1, N2, X1, M, Y1)
   integer, intent(IN) :: N1, N2, M
   double precision, intent(IN) :: X1(N1,N2)
   double precision, intent(OUT) :: Y1(N1)

   Y1(1:N1) = X1(1:N1, M)

end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine copy2LVecCX(N1, N2, X1, M, Y1)
   integer, intent(IN) :: N1, N2, M
   double complex, intent(IN) :: X1(N1,N2)
   double complex, intent(OUT):: Y1(N1)

   Y1(1:N1) = X1(1:N1, M)

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine copy3LVec(N1, N2, X1, M1, M2, Y1)
   integer, intent(IN) :: N1, N2, M1, M2
   double precision, intent(IN) :: X1(N1,N1,N2)
   double precision, intent(OUT) :: Y1(N1)

   Y1(1:N1) = X1(1:N1, M2, M1)

end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine copy3LVecCX(N1, N2, X1, M1, M2, Y1)
   integer, intent(IN) :: N1, N2, M1, M2
   double complex, intent(IN) :: X1(N1,N1,N2)
   double complex, intent(OUT):: Y1(N1)

   Y1(1:N1) = X1(1:N1, M2, M1)

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine copy4LHOSB(N0, N1, N2, X1, M1, M2, M3, Y1)
   integer, intent(IN) :: N0, N1, N2, M1, M2, M3
   double precision, intent(IN) :: X1(N0,N1,N1,N2)
   double precision, intent(OUT):: Y1(N0)

   Y1(1:N0) = X1(1:N0, M3, M2, M1)

end subroutine

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine copy4LHOSBCX(N0, N1, N2, X1, M1, M2, M3, Y1)
   integer, intent(IN) :: N0, N1, N2, M1, M2, M3
   double complex, intent(IN) :: X1(N0,N1,N1,N2)
   double complex, intent(OUT):: Y1(N0)

   Y1(1:N0) = X1(1:N0, M3, M2, M2)

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!*************************************************************
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function MA_dotProd(N, X1, Y1)
   integer, intent(IN)  :: N
   double precision, intent(IN) :: X1(N), Y1(N)

   MA_dotProd = dot_product(X1, Y1)

end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double complex function MA_dotProd_DX(N, X1, Y1)
   integer, intent(IN)  :: N
   double precision, intent(IN) ::  X1(N)
   double complex, intent(IN)   ::  Y1(N)

   double complex :: X(N)

   X(1:N)=X1(1:N)

   MA_dotProd_Dx = dot_product(X, Y1)

end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double complex function MA_dotProd_cx(N, X1, Y1)
   integer, intent(IN)  :: N
   double complex, intent(IN) :: X1(N), Y1(N)

   MA_dotProd_cx = dot_product(X1, Y1)

end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



!*************************************************************
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function HOSB_dotProd(B0, N, M, B, row, col, &
               Vi, Vj, HosbData)
   integer, intent(IN)  :: B0, N, M, B, row, col
   double precision, intent(IN) :: Vi(N), Vj(N)
   double precision, intent(IN) :: HosbData(N,M,M,B0)

   HOSB_dotProd = dot_product(Vi(1:N)*Vj(1:N),HosbData(1:N,row,col,B))

end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double complex function HOSB_dotProd_DX(B0, N, M, B, row, col, &
               Vi, Vj, HosbData)
   integer, intent(IN)  :: B0, N, M, B, row, col
   double precision, intent(IN) :: Vi(N), Vj(N)
   double complex, intent(IN)   :: HosbData(N,M,M,B0)

   double complex :: Vij(N)

   Vij(1:N) = DCMPLX(Vi(1:N)*Vj(1:N))

   HOSB_dotProd_DX = dot_product(Vij(1:N),HosbData(1:N,row,col,B))

end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double complex function HOSB_dotProd_CX(B0, N, M, B, row, col, &
               Vi, Vj, HosbData)
   integer, intent(IN)  :: B0, N, M, B, row, col
   double complex, intent(IN) :: Vi(N), Vj(N)
   double complex, intent(IN) :: HosbData(N,M,M,B0)

   HOSB_dotProd_CX = dot_product(Vi(1:N)*Vj(1:N),HosbData(1:N,row,col,B))

end function 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

