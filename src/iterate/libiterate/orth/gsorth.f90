!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             Grant-Schimit Orthogonization                  c
!c                  Just for one step                         c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

logical function GS_ORTH(N, WJ, M, VJ)
   implicit none
   integer, intent(IN) :: N, M
   double precision, intent(IN)  :: WJ(N)        
   double precision, intent(INOUT) :: VJ(N,M+1)   

!ccccccccccccccccc
   double precision :: dotwv
   integer :: i, M0

   M0 = M + 1
   VJ(1:N, M0) = WJ(1:N)   

   do I = 1, M
       dotwv = dot_product(WJ(1:N), VJ(1:N, I))
       VJ(1:N, M0) = VJ(1:N, M0) - dotwv * VJ(1:N, I)
   end do
  
   dotwv = DSQRT(dot_product(VJ(1:N, M0), VJ(1:N, M0)))

   if (dotwv == 0.0D0) then
      GS_ORTH = .false.
   else
      VJ(1:N, M0) = VJ(1:N, M0) / dotwv
      GS_ORTH = .true.
   end if

end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Modified Grant-Schimit Orthogonization           c
!c                  Just for one step                         c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

logical function MGS_ORTH(N, WJ, M, VJ)
   implicit none
   integer, intent(IN) :: N, M
   double precision, intent(IN)  :: WJ(N)        
   double precision, intent(INOUT) :: VJ(N,M+1)   

!ccccccccccccccccc   
   double precision :: dotwv, normw
   integer :: i, M0  

   M0 = M + 1
   VJ(1:N, M0) = WJ(1:N)   

   normw = DSQRT(dot_product(WJ(1:N), WJ(1:N)))
   MGS_ORTH = .false.

   if (normw == 0.0D0)   return

   do I = 1, M
       dotwv = dot_product(VJ(1:N, M0),  VJ(1:N,I))              
       VJ(1:N, M0) = VJ(1:N, M0) - dotwv * VJ(1:N, I)
   end do
  
   dotwv = DSQRT(dot_product(VJ(1:N, M0), VJ(1:N, M0)))

   if (dotwv /= 0.0D0) then      
      VJ(1:N, M0) = VJ(1:N, M0) / dotwv
      MGS_ORTH = .true.
   end if

end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             Grant-Schimit Orthogonization                  c
!c                 used for Lanczos                           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

logical function LAN_GS_ORTH(N, WJ, M, VJ, beta)
   implicit none
   integer, intent(IN) :: N, M
   double precision, intent(IN)  :: WJ(N)        
   double precision, intent(INOUT) :: VJ(N,M+1)  
   double precision, intent(OUT) :: beta

!ccccccccccccccccc
   double precision :: dotwv
   integer :: i, M0

   M0 = M + 1
   VJ(1:N, M0) = WJ(1:N)   

   do I = 1, M
       dotwv = dot_product(WJ(1:N), VJ(1:N, I))
       VJ(1:N, M0) = VJ(1:N, M0) - dotwv * VJ(1:N, I)
   end do
  
   dotwv = DSQRT(dot_product(VJ(1:N, M0), VJ(1:N, M0)))
   beta  = dotwv

   if (dotwv == 0.0D0) then
      LAN_GS_ORTH = .false.
   else
      VJ(1:N, M0) = VJ(1:N, M0) / dotwv
      LAN_GS_ORTH = .true.
   end if

end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                                            c
!c           Modified Grant-Schimit Orthogonization           c
!c                  Used for Lanczos algorithm                c
!c                                                            c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

logical function LAN_MGS_ORTH(N, WJ, M, VJ, beta)
   implicit none
   integer, intent(IN) :: N, M
   double precision, intent(IN)  :: WJ(N)        
   double precision, intent(INOUT) :: VJ(N,M+1)   
   double precision, intent(OUT)  :: beta

!ccccccccccccccccc   
   double precision :: dotwv, normw
   integer :: i, M0  

   M0 = M + 1
   VJ(1:N, M0) = WJ(1:N)   

   normw = DSQRT(dot_product(WJ(1:N), WJ(1:N)))
   if (normw == 0.0D0)     return

   do I = 1, M
       dotwv = dot_product(VJ(1:N, M0),  VJ(1:N,I))              
       VJ(1:N, M0) = VJ(1:N, M0) - dotwv * VJ(1:N, I)
   end do
  
   dotwv = DSQRT(dot_product(VJ(1:N, M0), VJ(1:N, M0)))
   beta = dotwv

   if (dotwv == 0.0D0) then  
      LAN_MGS_ORTH = .false.
   else    
      VJ(1:N, M0) = VJ(1:N, M0) / dotwv
      LAN_MGS_ORTH = .true.
   end if

end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             Grant-Schimit Orthogonization                  c
!c                  For the whole matrix                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

logical function GS_FULL_ORTH(nRow, nCol, Mat)
   implicit none
   integer, intent(IN) :: nRow, nCol
   double precision, intent(INOUT)  :: Mat(nRow, nCol)      
                         
!ccccccccccccccccc
   double precision :: tmpV (nRow), normV
   logical          :: GS_ORTH
   integer :: i

   normV = DSQRT(dot_product(Mat(1:nRow, 1), Mat(1:nRow, 1)))

   GS_FULL_ORTH = .false.
   if (normV == 0.0D0)     return

   Mat(1:nRow,1) = Mat(1:nRow, 1)/normV

   do I = 2, nCol
       tmpV(1:nRow)  = Mat(1:nRow, I)
       if ( .not. GS_ORTH(nRow, tmpV, I-1, Mat) )   return
   end do

   GS_FULL_ORTH = .true. 

end 
!****************************************************************


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Modified Grant-Schimit Orthogonization               c
!c                  For the whole matrix                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

logical function MGS_FULL_ORTH(nRow, nCol, Mat)
   implicit none
   integer, intent(IN) :: nRow, nCol
   double precision, intent(INOUT)  :: Mat(nRow, nCol)   
                         
!ccccccccccccccccc
   double precision :: tmpV(nRow), normV
   logical          :: MGS_ORTH
   
   integer :: i

   normV = DSQRT(dot_product(Mat(1:nRow, 1), Mat(1:nRow, 1)))

   MGS_FULL_ORTH = .false.
   if (normV == 0.0D0)       return

   Mat(1:nRow,1) = Mat(1:nRow, 1)/normV

   do I = 2, nCol
       tmpV(1:nRow)  = Mat(1:nRow, I)
       if ( .not. MGS_ORTH(nRow, tmpV, I-1, Mat) )  return
   end do

   MGS_FULL_ORTH = .true.

end 
!****************************************************************

