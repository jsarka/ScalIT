!
! Test for Pjm function
!

program testGaunt
    implicit none
    integer :: j1,j2,j3,m1,m2,m3, opt, nmax
    double precision :: gaunt, gaunt3j, gaunt3jComp, tmp1, tmp, threej
    double precision, allocatable :: lnn(:)
    j1 = 4 ;     m1 = 1 ;    j2 = 4
    m2 = 0 ;     j3 = 8 ;    m3 = 1

    DO
       nmax = j1+j2+j3+1
       ! print *, 'nmax:', nmax
       nmax = nmax + 1
       allocate (lnn(nmax))
       call lnFn(nmax, lnn)
       ! print *, 'n!:', lnn

       print *, 'j1=',j1, '  m1=',m1
       print *, 'j2=',j2, '  m2=',m2
       print *, 'j3=',j3, '  m3=',m3

       !print *, 'lnn:', lnn
       tmp = gaunt(j1,m1,j2,m2,j3,-m3)
       print *, 'Gaunt Test:(-m3)', tmp
       
       tmp1 = gaunt3jComp(j1,m1,j2, m2, j3, -m3, lnn)
       print *, 'Gaunt :', tmp1 

       tmp = gaunt3j(j1,m1,j2,m2,j3,-m3)
       print *, 'Gaung 3j:', tmp
        
       print *, 'Error:', tmp-tmp1

       deallocate(lnn)
    
       print *, 'Input optione: 0 for exit, other for continue'
       read *, opt
       if (opt==0) exit

       print *, 'Input J1, M1, j2, m2, j3, m3'
       read *, j1, m1, j2, m2,j3, m3 

    END DO

end
