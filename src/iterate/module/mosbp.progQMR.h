!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c             QMR Testing Subroutine                            c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine progQMR()
     double precision :: t1, t2, re
     double precision :: X1(pLen(sF)),B(plen(sF))
     double complex   :: X1CX(pLen(sF)),BCX(plen(sF))
     
     integer :: iter

     if (id==rootID) then
        print *
        print *, '*****************************************'
        print *, '*        Testing  QMR method            *'
        print *, '*****************************************'
        t1 = MPI_WTime()
     end if

     if (sOSB==TOSBW)  then
        if (.NOT. initMOSBW()) return
     end if

     if (sCX) then
         call randVec_cx(pLen(sF),BCX)
         iter = OSB_QMRCX(pLen(sF), BCX, X1CX, re)
     else
         call random_number(B)
         iter = OSB_QMR(pLen(sF), B, X1, re)
     end if

     if (sOSB==TOSBW)  call finalMOSBW()

     if (id==rootID) then
         t2 = MPI_WTime()
         print *
         print *, '   Preconditioner Method:', sOSB 
         print *, '   Iteration Num:',iter, '. Relative Error:', re
         print *
         print *, '   Total Time for QMR Testing:(sec)', t2-t1
         print *, '*****************************************************'
         print *
     end if

end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




