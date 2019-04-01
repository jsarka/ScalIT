!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function createWH()
    integer :: i, info, lwork
    double precision, allocatable :: WHWork(:)
    integer, parameter :: NB=6    

    createWH=.false.
    lwork = NB*hwLen    

    allocate(WHWORK(lwork), stat=info)
    if (info/=0) return 

    HWR(1:hwLen,1:hwLen)=HW(1:hwLen, 1:hwLen)
   
    select case (sPC)
    case (TA0) 
         call DSYTRF('U',hwLen,HWR,hwLen,WIPIV,WHWork,lwork,info)
    
    case (TA1) 
         do i = 1, hwLen
            HWR(i,i)=HW(i, i)-sOSBW%mE0
         end do
         call DSYTRF('U',hwLen,HWR,hwLen,WIPIV,WHWork,lwork,info)

    case default 
         do i = 1, hwLen
            HWR(i,i)=HW(i, i)-sOSBW%mE0
         end do

         call DSYTRF('U',hwLen,HWR,hwLen,WIPIV,WHWork,lwork,info)

    end select

    deallocate(WHWORK)

    createWH = (info==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function createWHCX()

    createWHCX = createWHCX_DX(.true.)   

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function createWHDX()

   createWHDX = createWHCX_DX(.false.)   

end function

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function createWHCX_DX(isCX)
    logical, intent(IN)  :: isCX

    integer :: i, info, lwork
    double complex, allocatable :: WHWorkCX(:)
    integer, parameter :: NB=6    

    createWHCX_DX=.false.
  
    if (hwLen <= 0) return

    lwork = NB * hwLen

    allocate(WHWORKCX(lwork), stat=info)
    if (info/=0) return       

    if (isCX) then
       HWRCX(1:hwLen,1:hwLen)=HWCX(1:hwLen, 1:hwLen)
    else
       HWRCX(1:hwLen,1:hwLen)=HW(1:hwLen, 1:hwLen)
    end if

    select case (sPC)
    case (TA0) 
         call ZHETRF('U',hwLen,HWRCX,hwLen,WIPIV,WHWorkCX,lwork,info)

    case (TA1) 
         do i = 1, hwLen
            HWRCX(i,i) = HWRCX(i,i)-sOSBW%mE0
         end do
         call ZHETRF('U',hwLen,HWRCX,hwLen,WIPIV,WHWorkCX,lwork,info)

    case (TAAP0) 
         do i = 1, hwLen
            HWRCX(i,i) = HWRCX(i,i)+DCMPLX(0.0D0,gAP(i))
         end do
         call ZGETRF(hwLen,hwLen,HWRCX,hwLen,WIPIV,info)

    case (TMAP0) 
         do i = 1, hwLen
            HWRCX(i,i) = HWRCX(i,i)+DCMPLX(0.0D0,-gAP(i))
         end do
         call ZGETRF(hwLen,hwLen,HWRCX,hwLen,WIPIV,info)

    case (TAAP) 
         do i = 1, hwLen
            HWRCX(i,i) = HWRCX(i,i)+DCMPLX(-sOSBW%mE0,gAP(i))
         end do
         call ZGETRF(hwLen,hwLen,HWRCX,hwLen,WIPIV,info)

    case (TAAPP) 
         do i = 1, hwLen
  	    HWRCX(i,i) = HWRCX(i,i)+DCMPLX(-sOSBW%mE0,gAPP(i))
         end do
         call ZGETRF(hwLen,hwLen,HWRCX,hwLen,WIPIV,info)

    case (TAAPR) 
         do i = 1, hwLen
            HWRCX(i,i) = HWRCX(i,i)+DCMPLX(-sOSBW%mE0,gAPR(i))
         end do
         call ZGETRF(hwLen,hwLen,HWRCX,hwLen,WIPIV,info)

    case (TMAP) 
         do i = 1, hwLen
            HWRCX(i,i) = HWRCX(i,i)+DCMPLX(-sOSBW%mE0,-gAP(i))
         end do
         call ZGETRF(hwLen,hwLen,HWRCX,hwLen,WIPIV,info)

    case (TMAPP) 
         do i = 1, hwLen
            HWRCX(i,i) = HWRCX(i,i)+DCMPLX(-sOSBW%mE0,-gAPP(i))
         end do
         call ZGETRF(hwLen,hwLen,HWRCX,hwLen,WIPIV,info)

    case (TMAPR) 
         do i = 1, hwLen
            HWRCX(i,i) = HWRCX(i,i)+DCMPLX(-sOSBW%mE0,-gAPR(i))
         end do
         call ZGETRF(hwLen,hwLen,HWRCX,hwLen,WIPIV,info)

    case default
         do i = 1, hwLen
            HWRCX(i,i) = HWRCX(i,i)-sOSBW%mE0
         end do
         call ZHETRF('U',hwLen,HWRCX,hwLen,WIPIV,WHWorkCX,lwork,info)
    end select

    deallocate(WHWORKCX)

    createWHCX_DX = (info==0)

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


