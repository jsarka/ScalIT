!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function initOSBW()
   integer :: info

   print *, 'original window width:', sOSBW%mDe,' size:',sosbw%mCnt
   call adjustWinSize(sOSBW%mCnt, sOSBW%mE0, sOSBW%mDE)
   hwLen = getOSBWSize(sOSBW%mE0, sOSBW%mDE)   
   print *, 'final window width:', sOSBW%mDe, '  window size:', sosbw%mCnt
   initOSBW = .false.

   if (hwlen==0) then
      sPC = sPC - 2;  return   ! use OSBD
   else
      if (.NOT.allocHW()) return
   end if

   call getOSBWIndex(sOSBW%mE0, sOSBW%mDE, hwLen, wind)

   if (sCX) then
      if (sHW>=0) then
          call getOSBWHCX(hwLen, wind, HWCX)
          initOSBW = .true.
          if (sHW>0) initOSBW = saveHWCX()
      else
          initOSBW = loadHWCX() 
      end if

      if (initOSBW)  initOSBW = createWHCX()
   else
      if (sHW>=0) then
          call getOSBWH(hwLen, wind, HW)
          initOSBW = .true.
          if (sHW>0) initOSBW = saveHW()
      else
          initOSBW = loadHW() 
      end if

      if (initOSBW) then
         if  (sJOB/=JOB_BOUND) then 
             initOSBW = createWHDX()
         else
             initOSBW = createWH()
         end if
      end if
   end if
 
   !print *, 'OSBW:len=',hwLen
   !print *,'HW:', HW
   !print *, 'HWR:', HWR

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine finalOSBW()
 
   call deallocHW()

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function getOSBWSize(E0, DE)
   double precision, intent(IN) :: E0, DE
   
   integer :: i

   getOSBWSize = 0
   do i=1, myLen
      if (ABS(Eig0(i)-E0) <= DE)           &
         getOSBWSize = getOSBWSize+1
   end do   
end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getOSBWIndex(E0, DE, N, ind)
   double precision, intent(IN) :: E0, DE
   integer (kind=8), intent(IN)  :: N
   integer, intent(OUT) :: ind(N) 
   
   integer :: i, cnt

   cnt = 1
   do i=1, myLen
      if (ABS(Eig0(i)-E0) <= DE) then
         ind(cnt) = i;         cnt = cnt + 1
         if (cnt > N) return
      end if
   end do   
end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getOSBWE0Ind(E0, DE, N, ind, Eout)
   double precision, intent(IN) :: E0, DE
   integer (kind=8), intent(IN)  :: N
   integer, intent(OUT) :: ind(N)
   double precision, intent(OUT) :: Eout(N) 
   
   integer :: i, cnt

   cnt = 1
   do i=1, myLen
      if (ABS(Eig0(i)-E0) <= DE) then
         ind(cnt) = i ; Eout(cnt)=Eig0(i)
         cnt = cnt + 1
         if (cnt > N) return
      end if
   end do   
end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine adjustWinSize(count, E0, de)
   integer, intent(IN)  :: count
   double precision, intent(IN)    :: E0
   double precision, intent(INOUT) :: de

   double precision :: de1
   integer :: i,cnt0

   de1=0.0D0;  cnt0 = 0;
   do i = 1,mylen 
      if (ABS(Eig0(i)-E0) <= de) cnt0 = cnt0 + 1
   end do

   do while (cnt0/=count)
      if (cnt0<count) then
         de1=de; de=2.0D0*de
      else
         de = 0.50D0*(de1+de)
      end if
  
      cnt0 = 0
      do i = 1,mylen 
         if (ABS(Eig0(i)-E0) <= de) cnt0 = cnt0 + 1
      end do

   end do

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function allocHW()
   integer :: info
  
   allocHW = .false.
   if (hwLen<=0) return
 
   allocate(WIPIV(hwLen), WIND(hwLen), stat=info)
   if (info /= 0) return

   allocate(HW(hwLen,hwLen),HWR(hwLen,hwLen),               &
           HWTMP(hwLen,hwLen),WX(hwLen),WY(hwLen),stat=info)  
   if (info /=0) return

   allocate(HWCX(hwLen,hwLen), HWRCX(hwLen,hwLen),          &
            HWTMPCX(hwLen,hwLen),WXCX(hwLen),WYCX(hwLen),   &
            stat=info) 
   if (info==0) allocHW = .true.

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine deallocHW()
    if (allocated(HW))     deallocate(HW)
    if (allocated(HWR))    deallocate(HWR)
    if (allocated(HWTMP))  deallocate(HWTMP)
    if (allocated(WX))     deallocate(WX)
    if (allocated(WY))     deallocate(WY)

    if (allocated(HWCX))   deallocate(HWCX)
    if (allocated(HWRCX))  deallocate(HWRCX)
    if (allocated(HWTMPCX))   deallocate(HWTMPCX)
    if (allocated(WXCX))   deallocate(WXCX)
    if (allocated(WYCX))   deallocate(WYCX)

    if (allocated(WIPIV))  deallocate(WIPIV)
    if (allocated(WIND))   deallocate(WIND)
end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

