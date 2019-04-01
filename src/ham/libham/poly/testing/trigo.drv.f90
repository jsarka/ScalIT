!
! Testing for Trigo, trigonom 
!
PROGRAM TEST_trigo
   implicit none

   integer :: Nmax, N, ang
   double precision, allocatable :: x0(:),x(:),cnx1(:,:),cnx2(:,:), snx1(:,:),snx2(:,:),tmp(:)
   integer :: I, OPT = 1 

   DO WHILE (OPT /= 0)
       PRINT *, 'Input N (# of x value), Nmax( max. of N*a ) '
       read  *, N, Nmax
       print *, ' Input input data type: 1=angle, other:cos(angle)'
       read *, ang
       if ((N < 1) .or. (Nmax < 2))   cycle

       allocate(x(N), tmp(N), x0(N))
       print *
       if (ang==1) then
          print *, 'Input theta values: data size=',N 
       else
          print *, 'Input cos(theta) values[-1.0, 1.0]: data size=',N 
       end if
       read *, x(1:N) 

       allocate(cnx1(N,Nmax), snx1(N,nmax), cnx2(N,Nmax+1),snx2(N,nmax+1))
       if (ang==1) then
          call Trigo1(Nmax, N, x, cnx1, snx1)
          call Trigonom1(Nmax, N, x, cnx2, snx2)
          x0(1:N) = x(1:N)
       else
          call Trigo2(Nmax, N, x, cnx1, snx1)
          call Trigonom2(Nmax, N, x, cnx2, snx2)
          x0(1:N) = acos(x(1:N))
       end if

       DO I = 1, Nmax    
          tmp(1:N)=cos(i*x0(1:N))  
          print *  
          print *, 'cos and sin: n=',i
          print *, 'D(cos)=',tmp(1:N)-cnx1(1:N,i)
          tmp(1:N)=sin(i*x0(1:N))
          print *, 'D(sin)=',tmp(1:N)-snx1(1:N,i)
       END DO

       print *
       print *, 'n=0'
       print *, 'cos:1.0', cnx2(1:N, 1)
       print *, 'sin:0.0', snx2(1:N, 1)
       DO I = 2, Nmax+1    
          print *
          print *, 'N=',i-1
          tmp(1:N)=cos((i-1)*x0(1:N))
          print *, 'D(cos)=',tmp(1:N)-cnx2(1:N,i)
          tmp(1:N)=sin((i-1)*x0(1:N))
          print *, 'D(sin)=',tmp(1:N)-snx2(1:N,i)
       END DO

       deallocate(x,x0, tmp,cnx1,cnx2, snx1, snx2) 
       PRINT *, 'INPUT OPTION for loop: O FOR EXIT'
       READ *, OPT
    END DO
END
