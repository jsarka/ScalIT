!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c           Eigen Value Testing Subroutine                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine progLanEig(eig)
    double precision, intent(OUT):: eig(sConv%mNum)

    integer :: iter
    double precision  :: t1, t2
    double precision  :: maxRES

    if (id==rootID) then
       t1 = MPI_WTime()
       print *
       print *, '********************************************'
       print *, '*     Eigen Values Using Lanczos method    *'
       print *, '********************************************'

       print *
       print *, 'Central Energy:', sConv%mE0
       print *,  '# of Interested Energy:', sConv%mNum
       print *, 'Final Convergence tolerence:', sConv%mTol
    end if

    iter = LANCZOS_CONV(sConv%mNum, Eig, maxRes)

    if (id==rootID) then
        print *,'Max Lanczos Error', maxres, ' Eigen Values:'
        print 1000, eig(1:sConv%mNum)
        t2 = MPI_WTime()
        print *, 'Time to get Lanczos Eigen values:(sec)',(t2-t1)
        print *, '************************************************'
    end if
    
 1000 format (E25.15,2x,E25.15,2x,E25.15)

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine progLanEig_CX(eig)
    double precision, intent(OUT) :: eig(sConv%mNum)

    integer :: iter
    double precision  :: t1, t2
    double precision  :: maxRES

    if (id==rootID) then
       t1 = MPI_WTime()
       print *
       print *, '*****************************************'
       print *, '*     Eigen Values Using PIST method    *'
       print *, '*****************************************'

       print *
       print *, 'Central Energy:', sConv%mE0
       print *,  '# of Interested Energy:', sConv%mNum 
       print *, 'Final Convergence tolerence:', sConv%mTol

    end if

    iter = LANCZOS_CONV(sConv%mNum, Eig, maxRes)

    if (id==rootID) then
       print *,'Max Lanczos Error', maxres, ' Eigen Values:'
       print 1000, eig(1:sConv%mNum)
       t2 = MPI_WTime()
       print *, 'Time to get LANCZOS Eigen values:(sec)',(t2-t1)
       print *, '************************************************'
    end if

 1000 format (E25.15,2x,E25.15,2x,E25.15)

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


