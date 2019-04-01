!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      Some utility subroutines for complex vector    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!******************************************************
double precision function getMax_CX(nType, N, vec)
     implicit none
     integer, intent(IN)  :: nType, N
     double complex, intent(in) :: vec(N)

     double precision :: opVal(N)
     double precision :: getMax     

     call getOpCX(nType, N, vec, opVal)
     getMax_CX = getMax(N, opVal)
     
end function

double precision function getMax0_CX(nType, N, vec, ind)
     implicit none
     integer, intent(IN)  :: nType, N
     double complex, intent(in) :: vec(N)
     integer, intent(out) :: ind

     double precision :: opVal(N)
     double precision :: getMax0     

     call getOpCX(nType, N, vec, opVal)
     getMax0_CX = getMax0(N, opVal, ind)
     
end function

double precision function getDiffMax_CX(nType, N, vec1, vec2)
     implicit none
     integer, intent(IN)  :: nType, N
     double complex, intent(in) :: vec1(N), vec2(N)

     double precision :: val1(N), val2(N)
     double precision :: getDiffMax

     call getOpCX(nType, N, vec1, val1)
     call getOpCX(nType, N, vec2, val2)

     getDiffMax_CX = getDiffMax(N, val1, val2)

end function
!*********************************************************

double precision function getDiffMax0_CX(nType, N, vec1, vec2, ind)
     implicit none
     integer, intent(IN)  :: nType, N
     double complex, intent(in) :: vec1(N), vec2(N)
     integer, intent(OUT) :: ind        

     double precision :: val1(N), val2(N)
     double precision :: getDiffMax0

     call getOpCX(nType, N, vec1, val1)
     call getOpCX(nType, N, vec2, val2)

     getDiffMax0_CX = getDiffMax0(N, val1, val2, ind)

end function

!*******************************
double precision function getSepMax_CX(nType, N, vec1, vec2)
     implicit none
     integer, intent(IN)  :: nType, N
     double complex, intent(in) :: vec1(N), vec2(N)

     double precision :: val1(N), val2(N)
     double precision :: getSepMax

     call getOpCX(nType, N, vec1, val1)
     call getOpCX(nType, N, vec2, val2)

     getSepMax_CX = getSepMax(N, val1, val2)

end function
!*********************************************************

double precision function getSepMax0_CX(nType, N, vec1, vec2, ind)
     implicit none
     integer, intent(IN)  :: nType, N
     double complex, intent(in) :: vec1(N), vec2(N)
     integer, intent(OUT) :: ind        

     double precision :: val1(N), val2(N)
     double precision :: getSepMax0

     call getOpCX(nType, N, vec1, val1)
     call getOpCX(nType, N, vec2, val2)

     getSepMax0_CX = getSepMax0(N, val1, val2, ind)

end function

!*******************************

double precision function getMin_CX(nType, N, vec)
     implicit none
     integer, intent(IN)  :: nType, N
     double complex, intent(in) :: vec(N)
 
     double precision :: opVal(N)
     double precision :: getMin

     call getOpCX(nType, N, vec, opVal)
     getMin_CX = getMin(N, opVal)
end function

double precision function getMin0_CX(nType, N, vec, ind)
     implicit none
     integer, intent(IN)  :: nType, N
     double complex, intent(in) :: vec(N)
     integer, intent(out) :: ind
 
     double precision :: opVal(N)
     double precision :: getMin0

     call getOpCX(nType, N, vec, opVal)
     getMin0_CX = getMin0(N, opVal, ind)
end function

double precision function getDiffMin_CX(nType, N, vec1, vec2)
     implicit none
     integer, intent(IN)  :: nType, N
     double complex, intent(in) :: vec1(N), vec2(N)

     double precision :: val1(N), val2(N)
     double precision :: getDiffMin

     call getOpCX(nType, N, vec1, val1)
     call getOpCX(nType, N, vec2, val2)

     getDiffMin_CX = getDiffMin(N, val1, val2)

end function

double precision function getDiffMin0_CX(nType, N, vec1, vec2, ind)
     implicit none
     integer, intent(IN)  :: nType, N
     double complex, intent(in) :: vec1(N), vec2(N)
     integer, intent(OUT) :: ind

     double precision :: val1(N), val2(N)
     double precision :: getDiffMin0

     call getOpCX(nType, N, vec1, val1)
     call getOpCX(nType, N, vec2, val2)

     getDiffMin0_CX = getDiffMin0(N, val1, val2, ind)

end function
!******************************************


subroutine AscReOrder_CX(ntype, N, vec)
    implicit none
    integer, intent(in) :: nType, N
    double complex, intent(INOUT)  :: vec(N)

!ccccccccccccccccccc
    integer        :: vecInd(N)
    double complex :: newVec(N)

    call ReOrder0_CX('A', nType, N, vec, vecInd, newVec)
    vec(1:N) = newVec(1:N)

end subroutine

subroutine DesReOrder_CX(ntype, N, vec)
    implicit none
    integer, intent(in) :: nType, N
    double complex, intent(INOUT)  :: vec(N)

!ccccccccccccccccccc
    integer        :: vecInd(N)
    double complex :: newVec(N)

    call ReOrder0_CX('D', nType, N, vec, vecInd, newVec)
    vec(1:N) = newVec(1:N)

end subroutine
!*********************************************

subroutine ReOrder_CX(order, ntype, N, vec)
    implicit none
    character*1, intent(in) :: order
    integer,intent(in)      :: nType, N
    double complex, intent(INOUT) :: vec(N)

    integer  :: vecInd(N)
    double complex :: newVec(N)  

    call ReOrder0_CX(order, ntype, N, vec, vecInd, newVec)

    vec(1:N) = newVec(1:N)
end subroutine


subroutine ReOrder0_CX(order, nType, N, vec, vecInd, newVec)
    implicit none
    character*1, intent(in) :: order
    integer,intent(in)      :: nType, N
    double complex, intent(IN)  :: vec(N)
    integer, intent(OUT)        :: vecInd(N)
    double complex, intent(OUT) :: newVec(N)

    double precision  :: opVec(N)

    call getOpCX(nType, N, vec, opVec)

    if ((order == 'D') .or.(order == 'd')) then          ! descend order
        call DesOrder(N, opVec, vecInd)       
    else
        call AscOrder(N, opVec, vecInd)
    end if

    newVec(1:N) = vec(vecInd(1:N))

end subroutine


!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Get the maximum one less than the given value   c
!cccccccccccccccccccccccccccccccccccccccccccccccccccc

double complex function getMaxUnder_CX(nType, N, vec, v0)
     implicit none
     integer, intent(IN)          :: nType, N
     double complex, intent(in)   :: vec(N)
     double precision, intent(in) :: v0     
     
     double precision :: opVec(N)
     integer :: ind

     call getOpCX(nType, N, vec, opVec)     
     call getMaxUnder0(N, opVec, v0, ind)

     getMaxUnder_CX = vec(ind)

end function

double complex function getMaxUnder0_CX(nType, N, vec, v0, ind)
     implicit none
     integer, intent(IN)          :: nType, N
     double complex, intent(in)   :: vec(N)
     double precision, intent(in) :: v0
     integer, intent(out)         :: ind
     
     double precision, dimension(N) :: opVec
  
     call getOpCX(nType, N, vec, opVec)     
     call getMaxUnder0(N, opVec, v0, ind)

     getMaxUnder0_CX = vec(ind)

end function

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      Get the Minimum greater than given values         c 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double complex function getMinAbove_CX(nType, N, vec, v0)
     implicit none
     integer, intent(IN)          :: nType, N
     double complex, intent(in)   :: vec(N)
     double precision, intent(in) :: v0

     integer          :: ind
     double precision :: opVec(N)
     
     call getOpCX(nType, N, vec, opVec)     
     call getMinAbove0(N, opVec, v0, ind)
     
     getMinAbove_CX = vec(ind)

end function

double complex function getMinAbove0_CX(nType, N, vec, v0, ind)
     implicit none
     integer, intent(IN)          :: nType, N
     double complex, intent(in)   :: vec(N)
     double precision, intent(in) :: v0
     integer, intent(out)         :: ind

     double precision :: opVec(N)
     
     call getOpCX(nType, N, vec, opVec)     
     call getMinAbove0(N, opVec, v0, ind)
     
     getMinAbove0_CX = vec(ind)

end function

!***********************************************************

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Get the index of values that are closest to E0            c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getWindow_CX(E0, nType, N, vec, M, newVec)
    implicit none
    double precision, intent(IN) :: E0
    integer, intent(IN)          :: nType, N, M
    double complex, intent(in)   :: vec(N)
    double complex, intent(out)  :: newVec(M)

    integer  :: vecInd(M)

    call getWindow0_CX(E0, nType, N, vec, M, vecInd, newVec)

end 


subroutine getWindow0_CX(E0, nType, N, vec, M, vecInd, newVec)
    implicit none
    double precision, intent(IN) :: E0
    integer, intent(IN)          :: nType, N, M
    double complex, intent(in)   :: vec(N)
    integer  , intent(out)       :: vecInd(M)
    double complex, intent(out)  :: newVec(M)

    double precision :: opVec(N),opnewVec(N)

!ccccccccccccccc    
  
    call getOpCX(nType, N, vec, opVec)   
    call getWindow0(E0, N, opVec, M, vecInd, opnewVec)

    newVec(1:M) = vec(vecInd(1:M))

end 
!***************************************************

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Test whether Convergence has met               c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function testConv_CX(ETOL, nType, N, vec1, vec2)
     implicit none
     double precision, intent(IN) :: ETOL
     integer, intent(IN)          :: nType, N
     double complex, intent(in)   :: vec1(N), vec2(N)

     double complex   :: err(N)
     double precision :: opErr(N)
     double precision :: maxErr, getMax
     integer          :: ind
     logical          :: testConv     

     err(1:N) = vec1(1:N)-vec2(1:N)
     call getOpCX(nType, N, err, opErr) 
     maxErr   = getMax(N, opErr, ind) 

     if (maxErr > ETOL) then
        testConv_CX = .false.
     else
        testConv_CX = .true.
     end if    

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      Get the value according to the operation of the complex       c
!c                           Type:                                    c
!c    1: real           2: |real|       3: imaginary                  c
!c    4: |imaginary|    5: |real| + |imaginary|                       c
!c    6: |real|^2 + |imaginary|^2       7: real+imaginary             c
!c    8: |real+imaginary|               9: real-imaginary             c
!c    10: imaginary-real               11: |imaginary-real|           c 
!c    12: |real|-|imaginary|           13: |imaginary|-|real|         c
!c    14: ||imaginary|-|real||         15: |real|^2-|imaginary|^2     c
!c    16: |imaginary|^2-|real|^2       17: ||real|^2-|imaginary|^2|   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getOpCX(nType, N, vec, op)
     implicit none
     integer, intent(in)           :: nType, N
     double complex, intent(in)    :: vec(N)
     double precision, intent(out) :: op(N) 
  
     select case (nType)   ! real
     case (1)
           op(1:N) = dble(vec(1:N))

     case (2)
           op(1:N) = abs(dble(vec(1:N)))

     case (3)
           op(1:N) = DIMAG(vec(1:N))

     case (4)
           op(1:N) = abs(DIMAG(vec(1:N)))

     case (5)
           op(1:N) = abs(dble(vec(1:N))) + abs(DIMAG(vec(1:N)))

     case (6)
           op(1:N) = dble(vec(1:N))**2+ DIMAG(vec(1:N))**2

     case (7)
           op(1:N) = dble(vec(1:N))+DIMAG(vec(1:N))

     case (8)
           op(1:N) = abs(dble(vec(1:N))+DIMAG(vec(1:N)))

     case (9)
           op(1:N) = dble(vec(1:N))-DIMAG(vec(1:N))

     case (10)
           op(1:N) = DIMAG(vec(1:N)) - dble(vec(1:N))

     case (11)
           op(1:N) = abs(dble(vec(1:N))-DIMAG(vec(1:N)))

     case (12)
           op(1:N) = abs(dble(vec(1:N)))-abs(DIMAG(vec(1:N)))

     case (13)
           op(1:N) = abs(DIMAG(vec(1:N)))-abs(dble(vec(1:N)))

     case (14)
           op(1:N) = abs(abs(dble(vec(1:N)))-abs(DIMAG(vec(1:N))))

     case (15)
           op(1:N) = dble(vec(1:N))**2-DIMAG(vec(1:N))**2

     case (16)
           op(1:N) = DIMAG(vec(1:N))**2 - dble(vec(1:N))**2

     case (17)
          op(1:N) = abs(dble(vec(1:N))**2-DIMAG(vec(1:N))**2)

     case default  ! CASE (0)   
           op(1:N) = dble(vec(1:N))

     end select

end subroutine
!***************************************************************
