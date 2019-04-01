!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Special case for CG coefficients         c
!c  C(j1,j2,j;m,(K-m),K) (0=<k=<kmax) will be      c
!c  calculate to calculate matrix elements in      c
!c  tetraatomic molecules. kmax=JTol               c
!ccccccccccccccccccccccccccccccccccccccccccccccccccc     
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    # of non-zero C(j1,j2,j,m,(K-m),K)           c
!c   where |m|<=j1, |k-m|<=j2, 0<=K<=min(j,kmax)   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function getCGmkSize(kmax, j1, j2, j)
   implicit none
   integer,intent(IN) :: kmax, j1, j2, j
    
   integer :: k0, k, m, km
   logical :: isTA

   getCGmkSize = 0
      
   if (.NOT. isTA(j1, j2, j))   return

   k0 = min(j, kmax)  

   do k = 0, k0
      do m = -j1, j1
         km = ABS(k-m)
         if (.NOT. (km>j2))  getCGmkSize = getCGmkSize + 1     
      end do   
    end do
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   (m,K) index for indth non-zero C(j1,j2,j,m,   c
!c   (K-m),K) and. mkInd(1)=m, mkInd(2)=K.         c
!c   where |m|<=j1, |k-m|<=j2, 0<=K<=min(j,kmax)   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getCGmkIndex(kmax, j1, j2, j, N, mkInd)
   implicit none
   integer, intent(IN)  :: kmax, j1, j2, j, N
   integer, intent(OUT) :: mkInd(2, N)
    
   integer :: i0, k0, km, k, m
   logical :: isTA
      
   if (.NOT. isTA(j1, j2, j))   return

   k0 = min(j, kmax)
   i0 = 1
   do k = 0, k0
      do m = -j1, j1
         km = ABS(k-m)
         if (.NOT. (km>j2)) then                   
            mkInd(1, i0) = m
            mkInd(2, i0) = k
            i0 = i0 + 1
            if (i0 > N) return         
         end if
      end do   
    end do
end

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   index for C(j1,j2,j, m,(K-m), K) in non-zero  c
!c   C(j1,j2,j,x,(X-x),x)). It returns 0/-1 when   c
!c   it is zeros, or not in the range.             c 
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function getCGmkPos(kmax, j1, j2, j, jm,jK)
   implicit none
   integer, intent(IN)  :: kmax, j1, j2, j, jm, jK
    
   integer :: i0, k0, k, m, km
   logical :: isTA

   getCGmkPos = -1
      
   if (.NOT. isTA(j1, j2, j))   return
   if ((jK<0) .OR. (ABS(jm)>j1))  return

   k0 = min(j, kmax)
  
   getCGmkPos = 0; i0 = 1
   do k = 0, k0
      do m = -j1, j1
         km = ABS(k-m)
         if (.NOT. (km>j2)) then          
            if ((k==jK) .AND. (m==jm)) then          
                getCGmkPos = i0;
                return;
            else
                i0 = i0 + 1
            end if        
         end if
      end do   
    end do
end


!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   (m,K) index for indth non-zero C(j1,j2,j,m,   c
!c   (K-m),K) and. mkInd(1)=m, mkInd(2)=K.         c
!c   where |m|<=j1, |k-m|<=j2, 0<=K<=min(j,kmax)   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function getCGmkInd(kmax, j1, j2, j, ind, jm, jk)
   implicit none
   integer, intent(IN)  :: kmax, j1, j2, j, ind
   integer, intent(OUT) :: jm, jk
    
   integer :: k0, k, m, km, i0
   logical :: isTA

   getCGmkInd = .FALSE.
      
   if (.NOT. isTA(j1, j2, j))   return

   k0 = min(j, kmax) 
   i0 = 1
   do k = 0, k0
      do m = -j1, j1
         km = ABS(k-m)
         if (.NOT. (km>j2)) then 
            if (i0 == ind) then            
                jm = m;   jk = k 
                getCGmkInd = .TRUE.
                return    
            else
                i0 = i0 + 1
            end if    
         end if
      end do   
    end do
end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Test whether (j1, j2, j3) forms a triangle       c
!c   i.e.: |j1-j2|<=j3<=|j1+j2|,|j2-j3|<=j1<=|j2+j3|  c
!c         |j3-j1|<=j2<=|j3+j1|. only is enough       c
!c   It returns true when (j1,j2,j3) forms a triangle c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc

logical function isTA(j1, j2, j3)
   implicit none
   integer, intent(IN) :: j1, j2, j3

   isTA = ((.NOT.(j1 < ABS(j2-j3))) .AND. (.NOT.( j1 > (j2+j3))))
   
end

logical function isNotTA(j1, j2, j3)
   implicit none
   integer, intent(IN) :: j1, j2, j3

   isNotTA = ((j1 < ABS(j2-j3)) .OR. ( j1 > (j2+j3)))
    
end
