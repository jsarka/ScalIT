!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!         Subroutine to calculate Gaunt coefficients:                   c
!             integral (Y(j1m1)*Y(j2m2)*Y(j3m3))                        c
!                                                                       c
! double gaunt(j1,m1, j2, m2, j3, m3, lnn). Do non-zero checking        c
! double gauntComp(j1,m1,j2,m2,j3,m3, lnn). Don't do non-zero checking  c
! lnn[j1+j2+j3+2], storing n!, n=0!,1!,2!,...,(j1+j2+j3+1)!             c
! logical isNonZeroGaunt(j1,m1,j2,m2,j3,m3) Return .TRUE. if Gaunt      c
! coeff is non-zero                                                     c
!                                                                       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

logical function isNonZeroGaunt(j1,m1,j2,m2,j3,m3)
    implicit none
    integer, intent(IN) :: j1,m1,j2,m2,j3,m3
    integer :: j123

    j123=j1+j2+j3 

    isNonZeroGaunt = ( ((m1+m2+m3)==0) .AND. (j123/2*2 == j123) .AND. &
                       (j1>=ABS(m1))   .AND. (j2>=ABS(m2))      .AND. &
                       (j3>=ABS(m3))   .AND.                          &
                       (j1<=(j2+j3))   .AND. (j1>=ABS(j2-j3)) .AND.   &
                       (j2<=(j1+j3))   .AND. (j2>=ABS(j1+j3)) .AND.   &
                       (j3<=(j1+j2))   .AND. (j3>=ABS(j1-j2))) 
end

!**********************************************************************
logical function isZeroGaunt(j1,m1,j2,m2,j3,m3)
    implicit none
    integer, intent(IN) :: j1,m1,j2,m2,j3,m3
    integer :: j123

    j123=j1+j2+j3 

    isZeroGaunt = ( ((m1+m2+m3)/=0) .OR. (j123/2*2 /= j123) .OR. &
                    (j1>=ABS(m1))   .OR. (j2>=ABS(m2))      .OR. &
                    (j3>=ABS(m3))   .OR.                         &
                    (j1<=(j2+j3))   .OR. (j1>=ABS(j2-j3))   .OR. &
                    (j2<=(j1+j3))   .OR. (j2>=ABS(j1+j3))   .OR. &
                    (j3<=(j1+j2))   .OR. (j3>=ABS(j1-j2))) 
end

!*************************************************************

double precision  function gaunt(j1,m1,j2,m2,j3,m3)
    implicit none
    integer, intent(IN) :: j1,m1,j2,m2,j3,m3

    double precision :: lnn(j1+j2+j3+2)
    logical :: isNonZeroGaunt
    double precision :: gauntComp

    ! check none-zero cases
    if (isNonZeroGaunt(j1,m1,j2,m2,j3,m3)) then
       call lnFn((j1+j2+j3+2), lnn)
       gaunt = gauntComp(j1,m1,j2,m2,j3,m3,lnn)
    else
       gaunt = 0.0D0
    end if
end 

double precision  function gaunt3j(j1,m1,j2,m2,j3,m3)
    implicit none
    integer, intent(IN) :: j1,m1,j2,m2,j3,m3

    double precision :: lnn(j1+j2+j3+2)
    logical :: isNonZeroGaunt
    double precision :: gaunt3jComp

    ! check none-zero cases
    if (isNonZeroGaunt(j1,m1,j2,m2,j3,m3)) then
       call lnFn((j1+j2+j3+2), lnn)
       gaunt3j = gaunt3jComp(j1,m1,j2,m2,j3,m3,lnn)
    else
       gaunt3j = 0.0D0
    end if
end


!*************************************************************
!  This version doesn't work due to the error in the paper   c
!*************************************************************
double precision  function gauntComp(j1,m1,j2,m2,j3,m3,lnn)
    implicit none
    integer, intent(IN) :: j1,m1,j2,m2,j3,m3
    double precision, intent(IN) :: lnn(j1+j2+j3+2)
    double precision, parameter  :: PI = 3.1415926535897932D0   

    integer :: L0, kmin, kmax, k
    integer :: k0, k1,k2,k3,k4,k5
    double precision :: fact1, fact

    L0 =(j1+j2+j3)/2
    k0 = j3+m1-m2+L0; k1 = j3-j1-m2; k2 = j3-j2+m1; 
    k3 = j1+j2-j3;    k4 = j1-m1;    k5 = j2+m2
    kmin = max(-k1, -k2);    kmin = max(kmin, 0)
    kmax = min(k3, k4);      kmax = min(kmax, k5)

    fact = 0.5D0*( lnn(j1+j2-j3+1) + lnn(j2+j3-j1+1) + lnn(j3+j1-j2+1)   &
                 - lnn(j1+j2+j3+2) + lnn(j1+m1+1)    + lnn(j1-m1+1)      &
                 + lnn(j2+m2+1)    + lnn(j2-m2+1)    + lnn(j3+m3+1)      &
                 + lnn(j3-m3+1)    )
    fact = fact + lnn(L0+1) - lnn(L0-j1+1) - lnn(L0-j2+1) - lnn(L0-j3+1)

    !print *, 'fact', fact, ' kmin:', kmin, ' kmax:', kmax   
    gauntComp = 0.0D0
    do k=kmin, kmax
       fact1 = fact -(  lnn(k+1)    + lnn(k+k1+1) + lnn(k+k2+1)       &
                      + lnn(k3-k+1) + lnn(k4-k+1) + lnn(k5-k+1) )
       if (k/2*2 == k) then
          gauntComp = gauntComp + exp(fact1)
       else
          gauntComp = gauntComp - exp(fact1)
       end if

    end do
    
    gauntComp = gauntComp * DSQRT((2*j1+1)*(2*j2+1)*(2*j3+1)/(4.0*PI))
    
    if (k0/2*2 /= k0)  gauntComp = - gauntComp

end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Using threej symbol to calculate Gaunt coefficient    c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision  function gaunt3jComp(j1,m1,j2,m2,j3,m3,lnn)
    implicit none
    integer, intent(IN) :: j1,m1,j2,m2,j3,m3
    double precision, intent(IN) :: lnn(j1+j2+j3+2)
 
    double precision, parameter  :: PI = 3.1415926535897932D0

    double precision :: threejComp

    gaunt3jComp = DSQRT((j1+0.5)*(j2+0.5)*(2*j3+1.0)/PI)    &
                  * threejComp(j1,0, j2,0, j3,0, lnn)       & 
                  * threejComp(j1,m1,j2,m2,j3,m3,lnn)
end

