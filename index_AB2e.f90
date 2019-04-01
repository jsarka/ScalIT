!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! (j,K)<->index transformation. used for triatomic
!  molecules calculation: AB2. r=rB-B  
!  For even parity(parity=.TRUE.): K=0,1,...,min(j,JTol)
!  For odd parity(parity=.false.): K=1,2,...,min(j,JTol)
!  j=[0, jmax]
!  For the permutation of B2 is constant, the function 
!  is also odd or even under permutation of B2. This version 
!  of wave function is snti-ymmetrry for B2 permutation, which 
!  means that j+K must be odd
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function get3Size(parity, JTol, jmax)
    implicit none
    logical, intent(IN) :: parity
    integer, intent(IN) :: JTol, jmax

    integer :: j, k, jk, kmin, kmax, jp, jt

    if (parity) then    ! even parity
       jp=0;
    else                ! odd parity
       jp=1
    end if
    jt=JTol+jp      ! total parity = (-1)^(p+JTol)

    if (jt/2*2==jt) then  ! even total parity
       kmin = 0; 
    else
       kmin = 1; 
    end if

!    jp=0
    get3Size = 0
    do j = 0, jmax
        kmax = min(JTol, j)
        do k = kmin, kmax
           jk = j + jp
           if (jk/2*2 == jk)  get3Size = get3Size + 1
        end do
    end do

end

subroutine get3Index(parity, JTol, jmax, N, jInd, kInd)
    implicit none
    logical, intent(IN) :: parity
    integer, intent(IN) :: JTol, jmax, N
    integer, intent(OUT):: jInd(N),kInd(N)

    integer :: j,k,jk,kmin,kmax,size, jp, jt
    
    if (parity) then
       jp=0;
    else
       jp=1;
    end if

    jt=JTol+jp
    if (jt/2*2==jt) then
        kmin = 0
    else
        kmin = 1
    end if
!    jp=0
    size = 0
    do j=0, jmax
       kmax = min(JTol, j)
       do k = kmin, kmax
          jk = j + jp
          if (jk/2*2 == jk) then
             size = size + 1
             if (size > N)  return
             jInd(size) = j; kInd(size)=k
           end if
       end do
     end do
end


