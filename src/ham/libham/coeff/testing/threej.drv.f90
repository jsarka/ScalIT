!
! Test for Pjm function
!

program test3j

    integer :: j1,j2,j3,m1,m2,m3, opt
    double precision :: threej, sixj, gaunt, tmp

    j1 = 4
    m1 = 1
    j2 = 4
    m2 = 0
    j3 = 8
    m3 =-1

    DO

       print *, 'j1=',j1, '  m1=',m1
       print *, 'j2=',j2, '  m2=',m2
       print *, 'j3=',j3, '  m3=',m3

       tmp = threej(j1,m1,j2,m2,j3,m3)
       print *, '3j symbol:', tmp

!       tmp = sixj(j1,m1,j2,m2,j3,m3)
!       print *, '6j symbol:', tmp

!       tmp = gaunt(j1,m1,j2,m2,j3,m3)
!       print *, 'Gaunt Coeff:', tmp
    
       print *, 'Input optione: 0 for exit, other for continue'
       read *, opt
       if (opt==0) exit

       print *, 'Input J1, M1, j2, m2, j3, m3'
       read *, j1, m1, j2, m2,j3, m3 

    END DO

end
