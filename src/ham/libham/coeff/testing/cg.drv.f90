!
! Test for Pjm function
!

program test3j
    implicit none
    integer :: j1,j2,j3,m1,m2,m3, opt, nmax
    double precision :: cg, tmp1, tmp, threej, cgComp, threejComp
    double precision, allocatable :: lnn(:)
    j1 = 4;    m1 = 1
    j2 = 4;    m2 = 0
    j3 = 8;    m3 = 1

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

       tmp  = threej(j1,m1,j2,m2,j3,-m3)
       tmp1 = threejComp(j1,m1,j2,m2,j3,-m3,lnn)
       print *, '3j symbol:(-m3)', tmp, tmp-tmp1

       print *
       print *, 'm:',m1, m2, m3       
       tmp = cg(j1,m1,j2, m2, j3, m3)
       tmp1 = cgComp(j1,m1,j2,m2,j3,m3, lnn)
       print *, 'CG symbol:', tmp, tmp-tmp1 
       
       print *
       tmp = threej(j1,m1,j2,m2,j3,-m3)      
       tmp = tmp*sqrt(2*j3+1.0D0)
       nmax = j1-j2-m3
       if (nmax/2*2 == nmax) then
          print *, 'Error:', tmp-tmp1
       else
          print *, 'Error:', tmp+tmp1
       end if    

       deallocate(lnn)
    
       print *, 'Input optione: 0 for exit, other for continue'
       read *, opt
       if (opt==0) exit

       print *, 'Input J1, M1, j2, m2, j3, m3'
       read *, j1, m1, j2, m2,j3, m3 

    END DO

end
