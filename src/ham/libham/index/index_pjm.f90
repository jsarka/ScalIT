!c    Have been checked using pjm.drv.f90
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Subroutines to deal with Pjm                      c
!c   The data are stored as:                           c
!c    m=0, j=0,1,2,...,jmax                            c
!c    m=1, j=1,2,...,jmax                              c
!c    ..................                               c
!c    m=jmax,j=jmax                                    c
!c So the total size of Pjm = (jmax+1)+jmax+...+1      c
!c     = (jmax+1)(jmax+2)/2                            c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function getPjmSize(jmax)
    implicit none
    integer, intent(IN) :: jmax
    
    getPjmSize = (jmax+1)*(jmax+2)/2
end 

!ccccccccccccccccccccccccccccccccccccccc
!c  (i)-> (jmInd(2,i)=[j(i), m(i)])    c
!ccccccccccccccccccccccccccccccccccccccc
subroutine getPjmIndex(jmax, jmInd)
    implicit none
    integer, intent(IN) :: jmax
    integer, intent(OUT):: jmInd(2, (jmax+1)*(jmax+2)/2)
  
    integer :: mi, mj, len

    len = 1
    do mi = 0, jmax
       do mj = mi, jmax
          jmInd(1, len ) = mj
          jmInd(2, len) = mi
          len = len + 1
       end do
    end do  
end

!
! (j, m) -> (index), 0=<j<=jmax, 0=<m<=j,
! otherwise, it returns -1
!
! The last Pjm of mth pjm is in the position (jmax+1)+...+(jmax+1-m)
!     = (m+1)*(2jmax+2-m)/2, so Pjm is in the position
!     of (m+1)*(2jmax+2-m)/2 - (j-m)
integer function getPjmPos(jmax, j, m)
    implicit none
    integer, intent(IN)  :: jmax, j, m

    if ((j>jmax) .or. (j<0) .or. (m>j) .or. (m<0)) then
        getPjmPos = -1
    else
        getPjmPos = (m+1)*(2*jmax+2-m)/2 - (jmax-j)
    end if
end

!
! (ind) -> (j, m)
! if ind > (jmax+1)(jmax+2)/2 or ind <=0, returns -1
! else return the (j, m) of indth in the Pjm order
!

subroutine getPjmInd(jmax, ind, j, m)
    implicit none
    integer, intent(IN)  :: jmax, ind
    integer, intent(OUT) :: j, m

    integer :: jmSize

    do m = 0, jmax
       jmSize = (m+1)*(2*jmax+2-m)/2;     
       if ((ind<jmSize) .OR. (ind==jmSize)) exit 
    end do 

    j = jmax - (jmSize - ind)

end

