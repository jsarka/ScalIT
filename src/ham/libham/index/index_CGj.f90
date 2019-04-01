!
! Index for Clebsch-Gordon coefficients
!

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Special case for CG coefficients         c
!c  C(j1,j2,j;m,(K-m),K) (0=<k=<kmax) will be      c
!c  calculate to calculate matrix elements in      c
!c  tetraatomic molecules. kmax=JTol               c
!ccccccccccccccccccccccccccccccccccccccccccccccccccc     
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    # of non-zero C(j1,j2,j,m,(K-m),K)           c
!c    due to j1, j2, j, where                      c
!c    0<=j1<=j1max, 0<=j2<=j2max,|j1-j2|<=j<=j1+j2 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function getCGjSize(j1max,j2max,jmax)
   implicit none
   integer,intent(IN) :: j1max, j2max, jmax

   integer :: i1,i2, j0min, j0max

   getCGjSize = 0

   do i1 = 0, j1max
      do i2 = 0, j2max
         j0min = abs(i1-i2);
         j0max = min((i1+i2), jmax)         
         getCGjSize = getCGjSize + (j0max-j0min+1) 
      end do      
   end do   
end


!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   index for C(j1,j2,j, m,(K-m), K) in non-zero  c
!c   C(j1,j2,j,x,(X-x),x)). It returns 0/-1 when   c
!c   it is zeros, or not in the range.             c 
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function getCGjPos(j1max,j2max,jmax, j1, j2, j)
   implicit none
   integer, intent(IN)  :: j1max, j2max, jmax, j1, j2, j
    
   integer :: i0, i1, i2, i3, j0max, j0min
   logical :: isTA

   getCGjPos = -1
      
   if (.NOT. isTA(j1, j2, j))       return
   if ( (j1 > j1max) .OR. (j1 <0) ) return
   if ( (j2 > j2max) .OR. (j2 <0) ) return
   if ( (j >  jmax)  .OR. (j  <0) ) return
  
   getCGjPos = 0; i0 = 1
   do i1 = 0, j1max
      do i2 = 0, j2max
         j0min = ABS(i1-i2) 
         j0max = min((i1+i2), jmax)
         do i3 = j0min, j0max            
            if ((i1==j1) .AND. (i2==j2) .AND. (i3==j)) then          
               getCGjPos = i0;
               return;
            else
               i0 = i0 + 1
            end if        
         end do
      end do   
    end do
end


!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   (m,K) index for indth non-zero C(j1,j2,j,m,   c
!c   (K-m),K) and. mkInd(1)=m, mkInd(2)=K.         c
!c   where |m|<=j1, |k-m|<=j2, 0<=K<=min(j,kmax)   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getCGjInd(j1max, j2max, jmax, ind, j1, j2, j)
   implicit none
   integer, intent(IN)  :: j1max, j2max, jmax, ind
   integer, intent(OUT) :: j1, j2, j
    
   integer :: i0, i1, i2, i3, j0min, j0max
     
   i0  = 1
   do i1 = 0, j1max
      do i2 = 0, j2max
         j0min = ABS(i1-i2)
         j0max = min((i1+i2), jmax)         
         do i3 = j0min, j0max
            if (i0==ind) then          
               j1=i1; j2=i2; j=i3; 
               return;
            else
               i0 = i0 + 1
            end if        
         end do
      end do   
    end do
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   (m,K) index for indth non-zero C(j1,j2,j,m,   c
!c   (K-m),K) and. mkInd(1)=m, mkInd(2)=K.         c
!c   where |m|<=j1, |k-m|<=j2, 0<=K<=min(j,kmax)   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getCGjIndex(j1max, j2max, jmax, N, jInd)
   implicit none
   integer, intent(IN)  :: j1max, j2max, jmax, N
   integer, intent(OUT) :: jInd(3, N)  

   integer :: i0, i1, i2, i3, j0min, j0max

   i0 = 1   
   do i1 = 0, j1max
      do i2 = 0, j2max
         j0min = ABS(i1-i2)
         j0max = min((i1+i2), jmax)       
         do i3 = j0min, j0max
            jInd(1, i0) = i1
            jind(2, i0) = i2
            jind(3, i0) = i3 
            i0 = i0 + 1
            if (i0 > N)  return                       
         end do
      end do   
    end do

end






