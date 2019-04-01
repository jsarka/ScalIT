
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       Calculate the threej symbols                        c
!c       ( j1/2  j2/2  j/2  )      where ji, mi are integers c
!c       ( m1/2  m2/2  m/2  )                                c
!c For the integer threej symbols, ji, mi are even integers  c
!c for the half-integer 3j symbols, ji, mi are odd integers  c
!c so the arguments are integer and twice the true value.    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double precision function twoThreeJ(j1,m1,j2,m2,j,m)
    implicit none
    integer, intent(IN) :: j1, m1, j2, m2, j, m   

    integer :: k, kmin, kmax, jtol
    integer :: jm1add, jm1Sub, jm2add, jm2Sub, jm3Add, jm3Sub
    double precision :: factor

    jtol = j1 + j2 + j

    twoThreeJ = 0.0D0
    ! check whether ji>=0 , ji>=|mi| and m1+m2+m3 = 0, and jtol is even, 
    if (((m1+m2+m3) /= 0) .OR. (jtol/2*2 /= jtol)          &
        .OR. (j1 < ABS(m1))  .OR. (j1 < 0)                 &
        .OR. (j2 < ABS(m2))  .OR. (j2 < 0)                 &
        .OR. (j  < ABS(m) )  .OR. (j < 0)  )
       return

   !  check whether ji, mi are both integer or both half-integer
    jm1Add = j1 + m1
    jm1Sub = j1 - m1
    jm2Add = j2 + m2
    jm2Sub = j2 - m2
    jm3Add = j3 + m3
    jm3Sub = j3 - m3

    if ((jm1Add/2*2 /= jm1Add) .OR. (jm1Sub/2*2 /= jm1Sub) .OR. &
        (jm2Add/2*2 /= jm2Add) .OR. (jm2Sub/2*2 /= jm2Sub) .OR. &
        (jm3Add/2*2 /= jm3Add) .OR. (jm3Sub/2*2 /= jm3Sub))  
        return

    jtol   = jtol/2       ! real (j1+j2+j3)
    jm1Add = jm1Add/2     ! real j1+m1
    jm1Sub = jm1Sub/2     ! real j1-m1
    jm2Add = jm2Add/2     ! real j2+m2
    jm2Sub = jm2Sub/2     ! real j2-m2
    jm3Add = jm3Add/2     ! real j3+m3
    jm3Sub = jm3Sub/2     ! real j3-m3



       factor = 0.0
       factor = binom(j1,(j1+j2-j)/2) / binom((j1+j2+j+2)/2,(j1+j2-j)/2)
       factor = factor * binom(j2,(j1+j2-j)/2) / binom(j1,(j1-m1)/2)
       factor = factor / binom(j2,(j2-m2)/2) / binom(j,(j-m)/2)
       factor = sqrt(factor)
       
       zmin = max(0,j2+(j1-m1)/2-(j1+j2+j)/2,j1+(j2+m2)/2-(j1+j2+j)/2)
       zmax = min((j1+j2-j)/2,(j1-m1)/2,(j2+m2)/2)
       
       sum=0.0
       do z = zmin,zmax
          par=1
          if(2*(z/2)-int(2*(z/2.0)) /= 0) par=-1
          
sum=sum+par*binom((j1+j2-j)/2,z)*binom((j1-j2+j)/2,(j1-m1)/2-z)*&
               binom((-j1+j2+j)/2,(j2+m2)/2-z)
       end do
       
       cleb = factor*sum
    end if
  end 

  recursive function binom(n,r) result(res)
    implicit none
    integer :: n,r
    real(rk) :: res
    real(rk) :: tmp
    if(n==r .or. r==0) then
       res = 1.0
    else if (r==1) then
       res = real(n,rk)
    else
       res = real(n,rk)/real(n-r,rk)*binom(n-1,r)
    end if
  end function binom

  recursive function factorial(n) result(res)
    implicit none
    integer :: n
    real(rk) :: res
    if (n==0 .or. n==1) then
       res=1.0
    else
       res=n*factorial(n-1)
    end if
  end function factorial


