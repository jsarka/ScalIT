!
! Index for Wigner 3j symbols or Clebsch-Gondon coefficients
!
! Storage order : http://www.siam.org/journals/sisc/25-4/42293.html
!

integer function get3jLSize(L)
   implicit none
   integer, intent(IN) :: L
   
   get3jLSize = L*(274+L*(225+L*(85+L*(L+15))))/120+1

end

integer function get3jSize(jmax)
   implicit none
   integer, intent(IN) :: jmax(3)
  
   integer :: L 
   L = max(jmax(1)+jmax(2)-jmax(3), jmax(2)+jmax(3)-jmax(1), &
           jmax(3)+jmax(1)-jmax(2))

   get3jSize = L*(274+L*(225+L*(85+L*(L+15))))/120+1
   
end 

subroutine get3jIndex( jmax, jSize, jInd)
   implicit none
   integer, intent(IN) :: jmax(3), jSize
   integer, intent(OUT) :: jInd(6, jSize)

   integer :: jm(3,2), rjm(3,2), L, X, T, B, S
   integer :: j1,j2,j3,m1,m2, j3min,j3max, ind
   logical :: reOrder3jm

   do j1 = 0, jmax(1)
      do j2 = 0, jmax(2)
         j3min = ABS(j1-j2)
         j3max = min(jmax(3), j1+j2)
         do j3 = j3min, j3max
            do m1 = -j1, j1
               do m2 = -j2, j2
                  jm(1,1)=j1; jm(2,1)=j2; jm(3,1)=j3
                  jm(1,2)=m1; jm(2,2)=m2; jm(3,2)=-(m1+m2)

                  call reOrder3jm(jm, rjm)  
                  L = rjm(1,1)-rjm(2,1)+rjm(3,1)
                  X = rjm(1,1)-rjm(1,2)
                  T = rjm(3,1)+rjm(3,2)
                  B = rjm(2,1)-rjm(2,2)
                  S = rjm(2,1)+rjm(3,1)-rjm(1,1)
                  ind = L*(24+L*(50+L*(35+L*(10+L))))/120 +         &
                        X*(6+X*(11+X*(6+X)))/24 + T*(2+T*(3+T))/6    &
                        +B*(B+1)/2 + S + 1
                   print *, 'Ind=', ind, L, X, T, B, S
                   print *, rjm
                   print *, jm
                   print *
                  if ((ind > 0) .AND. (ind < jSize)) then
                     if ((jInd(1,ind)==0) .AND. (jInd(2,ind)==0) .AND. &
                         (jInd(3,ind)==0) .AND. (jInd(4,ind)==0) .AND. &
                         (jInd(5,ind)==0) .AND. (jInd(6,ind)==0) )   then
                         jInd(1,ind)=j1; jInd(2,ind)=j2;jInd(3,ind)=j3
                         jInd(4,ind)=m1; jInd(5,ind)=m2;jInd(6,ind)=-(m1+m2) 
                     end if    
                 end if                
               end do
            end do
         end do
      end do
   end do
end 

integer function get3jPos(jInd)
   implicit none
   integer, intent(IN) :: jInd(3,2)

   integer :: rjm(3,2), L, X, T, B, S
   logical :: reOrder3jm

   call reOrder3jm(jInd, rjm)  
   L = rjm(1,1)-rjm(2,1)+rjm(3,1)
   X = rjm(1,1)-rjm(1,2)
   T = rjm(3,1)+rjm(3,2)
   B = rjm(2,1)-rjm(2,2)
   S = rjm(2,1)+rjm(3,1)-rjm(1,1)
   get3jPos = L*(24+L*(50+L*(35+L*(10+L))))/120           &
            + X*(6+X*(11+X*(6+X)))/24 + T*(2+T*(3+T))/6   &
            + B*(B+1)/2 + S + 1
end

!
! Reorder (j(1),j(2),j(3)) to (rj(1),rj(2),rj(3))
! so rj(1)>=rj(3)>=rj(2)
! it returns .TRUE. for the even permutations
! or returns .FALSE. for the odd permutations
!
logical function reOrder3j(j, rj)
   implicit none
   integer, intent(IN)  :: j(3)
   integer, intent(OUT) :: rj(3)

   if ( j(1) >= j (2) ) then
      if ( j(1) >= j(3) ) then   
          rj(1) = j(1)
          if (j(2) > j(3)) then  ! j(1)=max,j(3)=min
             rj(2) = j(3); rj(3)=j(2)
             reOrder3j = .FALSE.
          else                   ! j(1)=max,j(2)=min
             rj(2:3) = j(2:3)
             reOrder3j = .TRUE.
          end if 
      else                       ! j(3)=max,j(2)=min
          rj(1) = j(3); rj(2)=j(2)
          rj(3) = j(1); reOrder3j = .FALSE.
      end if
   else
     if (j(1)>j(3)) then     ! j(2)=max,j(3)=min
        rj(1)=j(2); rj(2)=j(3)
        rj(3)=j(1); reOrder3j = .TRUE.
     else
        rj(2)=j(1)          
        if (j(2)>j(3)) then  ! j(1)=min,j(2)=max
           rj(1)=j(2); rj(3)=j(3)
           reOrder3j = .FALSE.
        else
           rj(1)=j(3); rj(3)=j(2)
           reOrder3j = .TRUE.
        end if
     end if 
   end if
end

!
!  Reorder jm(1:3,1:2) to rjm(1:3,1:2
!  so rjm(1,1)>=rjm(3,1)>=rjm(2,1)
!    and rjm(2,2)>=0, and 
!    if (rjm(2,2)==0), rjm(3,2)>=0
!
!
logical function reOrder3jm(j, rj)
   implicit none
   integer, intent(IN)  :: j(3,2)
   integer, intent(OUT) :: rj(3,2)
   
   integer :: jSum   

   if ( j(1,1) >= j (2,1) ) then
      if ( j(1,1) >= j(3,1) ) then   
          rj(1,1:2) = j(1,1:2)
          if (j(2,1) > j(3,1)) then  ! j(1)=max,j(2)=min
             rj(2,1:2) = j(3,1:2)
             rj(3,1:2) = j(2,1:2)
             reOrder3jm = .FALSE.
          else                   ! j(1)=max,j(3)=min
             rj(2,1:2) = j(2,1:2)
             rj(3,1:2) = j(3,1:2)
             reOrder3jm = .TRUE.
          end if
      else                       ! j(3)=max,j(2)=min
          rj(1,1:2) = j(3,1:2); rj(2,1:2)  = j(2,1:2)
          rj(3,1:2) = j(1,1:2); reOrder3jm = .FALSE.
      end if
   else
     if (j(1,1)>=j(3,1)) then     ! j(2)=max,j(3)=min
        rj(1,1:2)=j(2,1:2); rj(2,1:2)  = j(3,1:2)
        rj(3,1:2)=j(1,1:2); reOrder3jm = .TRUE.
     else
        rj(2,1:2)=j(1,1:2)            
        if (j(2,1)>j(3,1)) then  ! j(1)=min,j(2)=max
           rj(1,1:2)=j(2,1:2); rj(3,1:2)=j(3,1:2)
           reOrder3jm = .FALSE.
        else
           rj(1,1:2)=j(3,1:2); rj(3,1:2)=j(2,1:2)
           reOrder3jm = .TRUE.
        end if
     end if 
   end if

   ! print *, 'rj:', rj
   if ( (rj(2,2) < 0) .OR. ((rj(2,2) == 0) .AND. (rj(3,2)<0)) ) then
      jsum = sum(j(1:3,1))
      rj(1:3,2) = - rj(1:3,2)
      if (jsum/2*2 /= jsum) reOrder3jm = (.NOT. reorder3jm)
   end if
end

