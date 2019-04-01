!cccccccccccccccccccccccccccccccccccccc
!c       Test for CG, 3j, Gaunt       c
!cccccccccccccccccccccccccccccccccccccc
program test_cgj

    use threejmod 
    use gauntmod

    implicit none

    interface getGaunt
       Module Procedure getGauntCoeff0
       Module Procedure getGauntCoeff1
       Module Procedure getGauntCoeff2
    end interface

    interface get3j
       Module Procedure get3jCoeff0
       Module Procedure get3jCoeff1
       Module Procedure get3jCoeff2
    end interface

    interface getCG
       Module Procedure getCGCoeff0
       Module Procedure getCGCoeff1
       Module Procedure getCGCoeff2
    end interface

    integer :: jmax(3), jm(3,2)
    integer :: j1,j2,j3,m1,m2,m3
    integer :: jmi, i, opt, nmax
    double precision :: cg, threej, gaunt, cgComp, threejComp, gauntComp
 
    double precision :: dj3(4), dj4(4), dcg(4), dga(4)
    logical :: isZero3j, isZeroCG, isZeroGaunt
    double precision, allocatable :: lnn(:)

    DO      
       print *, 'Input jmax for CG, 3j, Gaunt: 3 number'
       read *, jmax(1:3)
       nmax = sum(jmax(1:3))+2

       allocate (lnn(nmax))
       call lnFn(nmax, lnn)

       print *, 'Initialize Gaunt and 3j module'       
       if (initGa(jmax) .AND. init3j(jmax) ) then       


       print *, 'Do the calculation'  
       do j1=0,jmax(1)
          do j2=0,jmax(2)
             do j3 = 0, jmax(3)
                do m1 = -j1, j1
                   do m2 = -j2, j2
                      do m3 = -j3, j3
                         jm(1,1)=j1; jm(2,1)=j2; jm(3,1)=j3
                         jm(1,2)=m1; jm(2,2)=m2; jm(3,2)=-m3
                         print *
                         print *, j1,m1,j2,m2,j3,m3

                         dj3(1) = threej(j1,m1,j2,m2,j3,-m3)
                         if (isZero3j(j1,m1,j2,m2,j3,-m3)) then
                            dj3(2) = 0.0D0
                         else
                            dj3(2) = threejComp(j1,m1,j2,m2,j3,-m3,lnn)
                         end if

                         dj3(3) = get3j(jm)
                         dj3(4) = get3j(jm(1:3,1), jm(1:3,2))
                         dj3(5) = get3j(j1,m1,j2,m2,j3,-m3)
                         print *, dj3

                         dcg(1) = cg(j1,m1,j2,m2,j3,m3)
                         if (isZeroCG(j1,m1,j2,m2,j3,m3)) then
                            dcg(2) = 0.0D0
                         else
                            dcg(2) = cgComp(j1,m1,j2,m2,j3,m3,lnn)
                         end if
                         jm(3,2)= m3

                         dcg(3) = getcg(jm)
                         dcg(4) = getcg(jm(1:3,1), jm(1:3,2))
                         dcg(5) = getcg(j1,m1,j2,m2,j3,-m3)


                         dga(1) = gaunt(j1,m1,j2,m2,j3,m3)
                         if (isZeroGaunt(j1,m1,j2,m2,j3,m3)) then
                            dga(2) = 0.0D0
                         else
                            dga(2) = gauntComp(j1,m1,j2,m2,j3,m3,lnn)
                         end if

                         dga(3) = getGaunt(jm)
                         dga(4) = getGaunt(jm(1:3,1), jm(1:3,2))
                         dga(5) = getGaunt(j1,m1,j2,m2,j3,m3)

!                         dga(3) = getGauntCoeff0(jm)
!                         dga(4) = getGauntCoeff1(jm(1:3,1), jm(1:3,2))
!                         dga(5) = getGauntCoeff2(j1,m1,j2,m2,j3,m3)
                         
                         do i = 1, 4
                             dj4(i)=dsqrt(1.0D0+j3+j3)*dj3(i)
                         end do
                         jmi = j1-j2+m3
                         if (jmi/2*2 /= jmi) dj4(1:4) = - dj4(1:4)


                         print *, dj4
                         print *, dcg
                         print *, dga

                      end do
                   end do
                end do
             end do
          end do
       end do
       print *, 'Finish Calculation'
  
       call finalGA() ; call final3j()
!       call final3j() ; call finalGA()
       print *, 'deallocate lnn'
  
  
       end if
    deallocate(lnn)

       print *, 'Input optione: 0 for exit, other for continue'
       read *, opt
       if (opt==0) exit

    END DO

end
