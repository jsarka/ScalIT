!cccccccccccccccccccccccccccccccccccccccccccccc
!c               Get the maximum              c
!cccccccccccccccccccccccccccccccccccccccccccccc

double precision function getMax(N, vec)
     implicit none
     integer, intent(IN)  :: N
     double precision, intent(in) :: vec(N)

     integer :: I
     double precision    :: val     
     
     getMax = vec(1)   
     do I = 2, N         
          val = vec(I) 
          if (getMax < val ) then
               getMax = val               
          end if         
     end do

end function

double precision function getMax0(N, vec, ind)
     implicit none
     integer, intent(IN)  :: N
     double precision, intent(in) :: vec(N)
     integer, intent(out) :: ind

     integer :: I
     double precision    :: val     
     
     getMax0 = vec(1)
     ind    = 1
     do I = 2, N         
          val = vec(I) 
          if (getMax0 < val ) then
               getMax0 = val
               ind    = i
          end if         
     end do

end function

double precision function getAbsMax(N, vec)
     implicit none
     integer, intent(IN)  :: N
     double precision, intent(in) :: vec(N)

     integer :: I
     double precision    :: val     
     
     getAbsMax = abs(vec(1))   
     do I = 2, N         
          val = abs(vec(I))
          if (getAbsMax < val ) then
               getAbsMax = val               
          end if         
     end do

end function

double precision function getAbsMax0(N, vec, ind)
     implicit none
     integer, intent(IN)  :: N
     double precision, intent(in) :: vec(N)
     integer, intent(out) :: ind

     integer :: I
     double precision    :: val     
     
     getAbsMax0 = abs(vec(1))
     ind    = 1
     do I = 2, N         
          val = abs(vec(I))
          if (getAbsMax0 < val ) then
               getAbsMax0 = val
               ind    = i
          end if         
     end do

end function
!****************************************************************

double precision function getDiffMax(N, vec1, vec2)
     implicit none
     integer, intent(IN)  :: N
     double precision, intent(in) :: vec1(N), vec2(N)         

     integer :: I
     double precision    :: val     
     
     getDiffMax = abs(vec1(1)-vec2(1))     
     do I = 2, N         
          val = abs(vec1(I)-vec2(I))
          if (getDiffMax < val ) then
               getDiffMax = val               
          end if         
     end do

end function

double precision function getDiffMax0(N, vec1, vec2, ind)
     implicit none
     integer, intent(IN)  :: N
     double precision, intent(in) :: vec1(N), vec2(N) 
     integer, intent(OUT) :: ind

     integer :: I
     double precision    :: val     
     
     getDiffMax0 = abs(vec1(1)-vec2(1))
     ind  = 1   
     do I = 2, N         
          val = abs(vec1(I)-vec2(I))
          if (getDiffMax0 < val ) then
               getDiffMax0 = val
               ind    = i
          end if         
     end do

end function
!******************************************************

double precision function getSepMax(N1, vec1, N2, vec2)
     implicit none
     integer, intent(IN)  :: N1, N2
     double precision, intent(in) :: vec1(N1), vec2(N2)

     integer :: I, J
     double precision    :: val, minMax
 
     do I = 1, N1
          minMax = abs(vec1(I)-vec2(1))
          do J = 2, N2         
              val = abs(vec1(I)-vec2(J))
              if (minMax > val ) then
                  minMax = val               
              end if            
          end do 
          if (I == 1) then
              getSepMax = minMax
          else
              if (getSepMax < minMax) then
                 getSepMax = minMax
              end if    
          end if
     end do

end function

!****************************************************************
double precision function getSepMax0(N1, vec1, N2, vec2, ind)
     implicit none
     integer, intent(IN)  :: N1, N2
     double precision, intent(in) :: vec1(N1),vec2(N2)  
     integer, dimension(N1), intent(OUT) :: ind
     
     integer :: I, J
     double precision    :: val, minMax
 
     do I = 1, N1
          minMax = abs(vec1(I)-vec2(1))
          ind(I) = 1
          do J = 2, N2         
              val = abs(vec1(I)-vec2(J))
              if (minMax > val ) then
                  minMax = val      
                  ind(I) = J         
              end if            
          end do 
          if (I == 1) then
              getSepMax0 = minMax
          else
              if (getSepMax0 < minMax) then
                 getSepMax0 = minMax
              end if    
          end if
     end do
     
end function
!******************************************************


!cccccccccccccccccccccccccccccccccccccc
!c          Get the Minimum           c 
!cccccccccccccccccccccccccccccccccccccc
double precision function getMin(N, vec)
     implicit none
     integer, intent(IN)  ::  N
     double precision, intent(in) :: vec(N)

     integer :: I
     double precision :: val

     getMin = vec(1)     
     do I = 2, N
           val = vec(I)
           if (getMin > val ) then
               getMin = val             
           end if
     end do   
end function

double precision function getMin0(N, vec, ind)
     implicit none
     integer, intent(IN)  ::  N
     double precision, intent(in) :: vec(N)
     integer, intent(out) :: ind

     integer :: I
     double precision :: val

     getMin0 = vec(1)
     ind    = 1
     do I = 2, N
           val = vec(I)
           if (getMin0 > val ) then
               getMin0 = val
               ind    = i
           end if
     end do   
end function


double precision function getAbsMin(N, vec)
     implicit none
     integer, intent(IN)  ::  N
     double precision, intent(in) :: vec(N)

     integer :: I
     double precision :: val

     getAbsMin = abs(vec(1))    
     do I = 2, N
           val = abs(vec(I))
           if (getAbsMin > val ) then
               getAbsMin = val               
           end if
     end do   
end function

double precision function getAbsMin0(N, vec, ind)
     implicit none
     integer, intent(IN)  ::  N
     double precision, intent(in) :: vec(N)
     integer, intent(out) :: ind

     integer :: I
     double precision :: val

     getAbsMin0 = abs(vec(1))
     ind    = 1
     do I = 2, N
           val = abs(vec(I))
           if (getAbsMin0 > val ) then
               getAbsMin0 = val
               ind    = i
           end if
     end do   
end function


double precision function getDiffMin(N, vec1, vec2)
     implicit none
     integer, intent(IN)  :: N
     double precision, intent(in) :: vec1(N), vec2(N)

     integer :: I
     double precision    :: val     
     
     getDiffMin = abs(vec1(1)-vec2(1))     
     do I = 2, N         
          val = abs(vec1(I)-vec2(I))
          if (getDiffMin > val ) then
               getDiffMin = val               
          end if         
     end do

end function

double precision function getDiffMin0(N, vec1, vec2, ind)
     implicit none
     integer, intent(IN)  :: N
     double precision, intent(in) :: vec1(N), vec2(N)
     integer, intent(OUT) :: ind

     integer :: I
     double precision    :: val     
     
     getDiffMin0 = abs(vec1(1)-vec2(1))
     ind  = 1   
     do I = 2, N         
          val = abs(vec1(I)-vec2(I))
          if (getDiffMin0 > val ) then
               getDiffMin0 = val
               ind    = i
          end if         
     end do

end function
!******************************************************

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c          Get the Min and Max one time           c 
!ccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine getMinMax(N, vec,minVal, maxVal)
     implicit none
     integer, intent(IN)  ::  N
     double precision, intent(in)  :: vec(N)
     double precision, intent(out) :: minVal, maxVal     

     integer :: I
     double precision :: val

     minVal = vec(1)
     maxVal = minVal  
     do I = 2, N
           val = vec(I)
           if (minVal > val ) then
               minVal = val               
               cycle
           end if
           if (maxVal < val ) then
               maxVal = val              
           end if
     end do   
end subroutine

subroutine getMinMax0(N, vec,minVal, minInd, maxVal, maxInd)
     implicit none
     integer, intent(IN)  ::  N
     double precision, intent(in)  :: vec(N)
     double precision, intent(out) :: minVal, maxVal
     integer, intent(out)   :: minInd, maxInd

     integer :: I
     double precision :: val

     minVal = vec(1)
     maxVal = minVal
     minInd = 1
     maxInd = 1
     do I = 2, N
           val = vec(I)
           if (minVal > val ) then
               minVal = val
               minInd    = i
               cycle
           end if
           if (maxVal < val ) then
               maxVal = val
               maxInd = i
           end if
     end do   
end subroutine


subroutine getAbsMinMax(N, vec,minVal, maxVal)
     implicit none
     integer, intent(IN)  ::  N
     double precision, intent(in)  :: vec(N)
     double precision, intent(out) :: minVal, maxVal     

     integer :: I
     double precision :: val

     minVal = abs(vec(1))
     maxVal = minVal 
     do I = 2, N
           val = abs(vec(I))
           if (minVal > val ) then
               minVal = val             
               cycle
           end if
           if (maxVal < val ) then
               maxVal = val            
           end if
     end do   
end subroutine

subroutine getAbsMinMax0(N, vec,minVal, minInd, maxVal, maxInd)
     implicit none
     integer, intent(IN)  ::  N
     double precision, intent(in)  :: vec(N)
     double precision, intent(out) :: minVal, maxVal
     integer, intent(out)   :: minInd, maxInd

     integer :: I
     double precision :: val

     minVal = abs(vec(1))
     maxVal = minVal
     minInd = 1
     maxInd = 1
     do I = 2, N
           val = abs(vec(I))
           if (minVal > val ) then
               minVal = val
               minInd    = i
               cycle
           end if
           if (maxVal < val ) then
               maxVal = val
               maxInd = i
           end if
     end do   
end subroutine

!******************************************************

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c               Reorder the vector                c 
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine  Reorder (order, N, vec)
     implicit none
     character*1, intent(IN) :: order    ! 'A'/'a'/'D'/'d'
     integer, intent(IN)     :: N
     double precision, intent(INOUT) :: vec(N)

     integer          :: vecInd(N)
     double precision :: newVec(N)

     call Reorder0(order, N, vec, vecInd, newVec)

     vec(1:N) = newVec(1:N)

end subroutine

!*****************************************************
subroutine  Reorder0 (order, N, vec, vecInd, newVec)
     implicit none
     character*1, intent(IN) :: order ! 'A'/'a'/'D'/'d'
     integer, intent(IN)     :: N
     double precision, intent(IN)  :: vec(N)
     integer, intent(OUT)          :: vecInd(N)
     double precision, intent(OUT) :: newVec(N)

     integer :: I
   
     if (order=='D' .or. order=='d') then
          call DESORDER(N, vec, vecInd)
     else
          call ASCORDER(N, Vec, vecInd)
     end if  

     do I = 1, N        
         newVec(I) = vec(vecInd(I))
     end do

end subroutine

!******************************************************

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c          Reorder the vector in Asc order        c 
!ccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine AscOrder( N, vec, vecInd)     ! Ascend order
     implicit none     
     integer, intent(IN)      :: N
     double precision, intent(IN) :: vec(N)
     integer, intent(out)     :: vecInd(N)
 
     integer :: i,j, index
     logical, dimension(N) :: label
     double precision      :: val0, val1

     vecInd(1:N) = 0
     label(1:N)  = .true.

     do I = 1, N
        index = 0
        do J = 1, N
           if (label(J)) then          ! need to deal with
               if (index == 0) then    ! initial values
                   index = J
                   val0  = vec(J)
               else
                   val1 = vec(J)
                   if (val0 > val1) then
                       val0 = val1
                       index = J
                   end if    ! finish val0 >= val1
               end if   ! finish index /= 0
           end if       ! finish label
        end do          ! finish one loop
        vecInd(I) = index
        label(index) = .false.
     end do
       
end subroutine


!******************************************************

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c          Reorder the vector in Des order        c 
!ccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine DesOrder(N, vec, vecInd)    ! Descend Order
     implicit none     
     integer, intent(IN)      :: N
     double precision, intent(IN) :: vec(N)
     integer, intent(out)     :: vecInd(N)
 
     integer :: i,j, index
     logical :: label(N)
     double precision      :: val0, val1

     vecInd(1:N) = 0
     label(1:N)  = .true.

     do I = 1, N
        index = 0
        do J = 1, N
           if (label(J)) then          ! need to deal with
               if (index == 0) then    ! initial values
                   index = J
                   val0  = vec(J)
               else
                   val1 = vec(J)
                   if (val0 < val1) then
                       val0 = val1
                       index = J
                   end if    ! finish val0 >= val1
               end if   ! finish index /= 0
           end if       ! finish label
        end do          ! finish one loop
        vecInd(I) = index
        label(index) = .false.
     end do

end subroutine
!**********************************************************


!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Get the maximum one less than the given value   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function getMaxUnder(N, vec, v0)
     implicit none
     integer, intent(IN)          ::  N
     double precision, intent(in) :: vec(N)
     double precision, intent(in) :: v0     

     integer :: I
     double precision :: val, minVal, maxVal    
       
     call getMinMax(N, vec, minVal, maxVal)

     if ((V0 <= maxVal) .and. (v0 >= minVal)) then 
        
         do I = 1, N
           val = vec(I)
           if (( val > minval ) .and. (val < V0) ) then
               minVal = val              
           end if
         end do
     end if
     
     getMaxUnder = minVal

end function

double precision function getMaxUnder0(N, vec, v0, ind)
     implicit none
     integer, intent(IN)          ::  N
     double precision, intent(in) :: vec(N)
     double precision, intent(in) :: v0
     integer, intent(out)         :: ind

     integer :: I, minInd, maxInd
     double precision :: val, minVal, maxVal    
       
     call getMinMax0(N, vec, minVal, minInd, maxVal, maxInd)

     ind = 0     

     if ((V0 <= maxVal) .and. (v0 >= minVal)) then 
         ind = minInd
         do I = 1, N
           val = vec(I)
           if (( val > minval ) .and. (val < V0) ) then
               minVal = val
               ind    = i
           end if
         end do
     end if
     
     getMaxUnder0 = minVal

end function
!***********************************************************

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      Get the Minimum greater than given values         c 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function getMinAbove( N, vec, v0)
     implicit none
     integer, intent(IN)          :: N
     double precision, intent(in) :: vec(N)
     double precision, intent(in) :: v0    

     integer :: I
     double precision :: val, minVal, maxVal     
     
     call getMinMax(N, vec, minVal,maxVal)
     
     if ((V0 <= maxVal) .and. (v0 >= minVal)) then           
         do I = 1, N
           val = vec(I)
           if (( val < maxval ) .and. (val > V0) ) then
               maxVal = val               
           end if
         end do
     end if
     
     getMinAbove = maxVal

end function

double precision function getMinAbove0( N, vec, v0, ind)
     implicit none
     integer, intent(IN)          :: N
     double precision, intent(in) :: vec(N)
     double precision, intent(in) :: v0
     integer, intent(out)         :: ind

     integer :: I, minInd, maxInd
     double precision :: val, minVal, maxVal     
     
     call getMinMax0(N, vec, minVal, minInd, maxVal, maxInd)

     ind = 0

     if ((V0 <= maxVal) .and. (v0 >= minVal)) then  
         ind = maxInd    
         do I = 1, N
           val = vec(I)
           if (( val < maxval ) .and. (val > V0) ) then
               maxVal = val
               ind    = i
           end if
         end do
     end if
     
     getMinAbove0 = maxVal

end function
!**********************************************************

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Get the index of values that are closest to E0            c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getWindow(E0, N, vec, M, newVec)
    implicit none
    double precision, intent(IN) :: E0
    integer, intent(IN)          :: N, M
    double precision, intent(in) :: vec(N)
    double precision, intent(out):: newVec(M)

    integer,dimension(M)  :: vecInd

    call getWindow0(E0, N, vec, M, vecInd, newVec)
end 

!**********************************************************
subroutine getWindow0(E0, N, vec, M, vecInd, newVec)
    implicit none
    double precision, intent(IN) :: E0
    integer, intent(IN)          :: N, M
    double precision, intent(in) :: vec(N)
    integer  , intent(out)       :: vecInd(M)
    double precision, intent(out) :: newVec(M)

    double precision :: absVec(N)
    logical          :: label(N)
    double precision :: val0, val1
    integer :: I, J, index  

!ccccccccccccccc    
    if (N<M) then
       print *, 'Warning, N < M in calling getWindow0'
       newVec(1:N)=vec(1:N); newVec(N:M)=vec(N)
       vecInd(N:M)=N
       do i = 1, N
          vecInd(i)=i
       end do
       return
    end if
 
    absVec(1:N) = abs(vec(1:N)-E0)
    vecInd(1:M) = 0;    label(1:N)  = .true.

     do I = 1, M
        index = 0
        do J = 1, N
           if (label(J)) then          ! need to deal with
               if (index == 0) then    ! initial values
                   index = J
                   val0  = absVec(J)
               else
                   val1 = absVec(J)
                   if (val0 > val1) then
                       val0 = val1
                       index = J
                   end if    ! finish val0 >= val1
               end if   ! finish index /= 0
           end if       ! finish label
        end do          ! finish one loop
        vecInd(I) = index
        newVec(I) = vec(index)
        label(index) = .false.
     end do
end 
!***************************************************
!**********************************************************
subroutine getWindow1(E0, N, vec, M, vecInd)
    implicit none
    double precision, intent(IN) :: E0
    integer, intent(IN)          :: N, M
    double precision, intent(in) :: vec(N)
    integer  , intent(out)  :: vecInd(M)

    double precision :: absVec(N)
    logical          :: label(N)
    double precision :: val0, val1
    integer :: I, J, index  

!ccccccccccccccc    
    if (N<M) then
       print *, 'Warning, N < M in calling getWindow1'
       vecInd(N:M)=N
       do i = 1, N
          vecInd(i)=i
       end do
       return
    end if

    absVec(1:N) = abs(vec(1:N)-E0)
    vecInd(1:M) = 0;    label(1:N)  = .true.

     do I = 1, M
        index = 0
        do J = 1, N
           if (label(J)) then          ! need to deal with
               if (index == 0) then    ! initial values
                   index = J
                   val0  = absVec(J)
               else
                   val1 = absVec(J)
                   if (val0 > val1) then
                       val0 = val1
                       index = J
                   end if    ! finish val0 >= val1
               end if   ! finish index /= 0
           end if       ! finish label
        end do          ! finish one loop
        vecInd(I) = index
        label(index) = .false.
     end do

end 



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Test whether Convergence has met               c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function testConv(ETOL, N, vec1, vec2)
     implicit none
     double precision, intent(IN) :: ETOL
     integer, intent(IN)          :: N
     double precision, intent(in) :: vec1(N), vec2(N)

     double precision, dimension(N)  :: err
     double precision :: getMax0, maxErr 
     integer          :: ind   

     err(1:N) = abs(vec1(1:N)-vec2(1:N))
     maxErr  = getMax0(N, err, ind) 

     if (maxErr > ETOL) then
        testConv = .false.
     else
        testConv = .true.
     end if

end function

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Get the index of values that are closest to E0            c
!c    N >= M, Doesn't work when some values are the same         c
!c    If the values are the same, only the first one is chosen   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getWindow_EX(N, vec,E0, M, ind, newVec)
    implicit none
    double precision, intent(IN) :: E0
    integer, intent(IN)          :: N, M
    double precision, intent(in) :: vec(N)
    integer  , intent(out)       :: ind(M), newVec(M)

    double precision :: absVec(N) 
    integer :: I
    double precision :: val
    double precision :: getMin0, getMinAbove0

    absVec(1:N) = abs(vec(1:N)-E0)
    val = getMin0(N, absVec, ind(1))     
    newVec(1) = vec(ind(1))
   
    do I = 2, M               
        val = getMinAbove0(N, absVec, val, ind(I))
        newVec(I) = vec(ind(I))
    end do      

end subroutine

!********************************************************
