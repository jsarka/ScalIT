
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine adjustWinSize(count, E0, de)
   integer, intent(IN)  :: count
   double precision, intent(IN)    :: E0
   double precision, intent(INOUT) :: de

   double precision :: de1
   integer :: i,cnt0, allCount, ierr

   de1=0.0D0;  cnt0 = 0;
   do i = 1,pLen(1)
      if (ABS(Eig0(i)-E0) <= de) cnt0 = cnt0 + 1
   end do

   call MPI_ALLREDUCE(cnt0, allCount, 1, MPI_INTEGER, MPI_SUM,   &
                         MPI_COMM_WORLD,IERR)   

   do while (allCount/=count)
      if (allCount<count) then
         de1=de; de=2.0D0*de
      else
         de = 0.50D0*(de1+de)
      end if
  
      cnt0 = 0
      do i = 1,pLen(1) 
         if (ABS(Eig0(i)-E0) <= de) cnt0 = cnt0 + 1
      end do

      call MPI_ALLREDUCE(cnt0, allCount, 1, MPI_INTEGER, MPI_SUM,   &
                         MPI_COMM_WORLD,IERR) 
   end do

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getWGlobalSize(E0, DE, wLcnt, wGCnt)
      double precision, intent(IN) :: E0, DE
      integer, intent(OUT) :: wLCnt,wGCnt

      integer :: i, ierr
      double precision :: Emin, EMax

      Emin = Min(E0-DE,E0+DE);  Emax=Max(E0-DE,E0+DE)

      wLCnt = 0; wGCnt = 0

      do i = 1, plen(1)
         if ((Eig0(i)>=Emin) .AND. (Eig0(i)<=Emax))   &
             wLCnt = wLCnt + 1
      end do  

      call MPI_ALLREDUCE(wLCnt, wGCnt, 1, MPI_INTEGER, MPI_SUM,    &
                         MPI_COMM_WORLD,IERR)

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getWGlobalIndex(E0, DE, cnt1, cnt2, glbcnt, sumCnt, locInd, &
                           glbInd, nodInd )
   double precision, intent(IN) :: E0, DE
   integer, intent(IN)  :: cnt1, cnt2
   integer, intent(OUT) :: glbCnt(nNodes), sumCnt(nNodes)
   integer, intent(OUT) :: locInd(cnt1),glbInd(cnt2),nodInd(cnt2)

   integer :: i, len, cnt, ierr, locNind(cnt2)

   double precision :: Emin, Emax

   Emin = Min(E0-DE,E0+DE);  Emax=Max(E0+DE,E0-DE)

   len = 1; cnt=0
   do i = 1, pLen(1)
       if ((Eig0(i)>=Emin) .AND. (Eig0(i)<=Emax))  then
           cnt=cnt+1
           if (len<=cnt1) then 
  	       locind(len)=i;
               locNind(len)=id
               len=len+1               
           end if
       end if
   end do 
  
   call MPI_ALLGATHER(cnt, 1, MPI_INTEGER, glbCnt, 1, MPI_INTEGER,   &
                         MPI_COMM_WORLD, IERR)
   sumCnt(1) = 0       
   do I = 1, nNodes-1
       sumCnt(I+1) = sumCnt(I) + glbCnt(I)
   end do
          
   call MPI_ALLGATHERV(locInd, cnt, MPI_INTEGER, glbInd, glbCnt, sumCnt,  &
                          MPI_INTEGER, MPI_COMM_WORLD, IERR) 
 
   call MPI_ALLGATHERV(locNind, cnt, MPI_INTEGER, nodInd, glbCnt, sumCnt, &
                          MPI_INTEGER, MPI_COMM_WORLD, IERR)  

end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calGlobalIndex(nSize, gpInd, ggInd, giInd)
   integer, intent(IN) :: nSize, gpInd(hwLen)
   integer(kind=MPI_OFFSET_KIND), intent(OUT) :: ggInd(hwLen)
   integer, intent(OUT) :: gIInd(hwLen)

   integer :: i
   integer(kind=MPI_OFFSET_KIND) :: offind(nNodes), pnum, qnum   

   gIInd(1:hwLen) = 0
   if (phwLen>0)  giInd(sCnt(id+1)+1:sCnt(id+1)+phwLen)=pInd(1:phwLen)

   pnum = myconf%gBlk(1)/nNodes
   qnum = myconf%gBlk(1) - pnum * nNodes
   offind(1) = 0

   do i = 1, qnum
      offind(i+1) = offind(i)+(pnum+1)*sN(1)
   end do

   do i = qnum+1, nNodes-1  
      offind(i+1) = offind(i) + pnum*sN(1)
   end do
   
   do i = 1, hwlen
     ggInd(i) = gInd(i)+offind(nodInd(i)+1)
   end do

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 integer function getWinSize(E0, DE)
    double precision, intent(IN) :: E0,DE 
    integer :: len    

    call getWGlobalSize((E0-DE),(E0+DE),len,getWinSize)
         
 end function 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
integer function getWLocalSize(E0, DE)
      double precision, intent(IN) :: E0, DE

      integer :: i
      double precision :: Emin, EMax

      Emin = Min(E0-DE,E0+DE);  Emax=Max(E0-DE,E0+DE)

      getWLocalSize = 0

      do i = 1, plen(1)
         if ((Eig0(i)>=Emin) .AND. (Eig0(i)<=Emax))  &
              getWLocalSize = getWLocalSize + 1
      end do   
      
end function
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getWLocalIndex(E0, DE, cnt, ind)
      double precision, intent(IN) :: E0, DE
      integer, intent(IN)  :: cnt
      integer, intent(OUT) :: ind(cnt)

      integer :: i, len
      double precision :: Emin, EMax

      Emin = Min(E0-DE,E0+DE);  Emax=Max(E0-DE,E0+DE)
      len = 1
      do i = 1, pLen(1)
         if ((Eig0(i)>=Emin) .AND. (Eig0(i)<=Emax))  then
             ind(len)=i; len=len+1
             if (len>cnt) return
         end if
      end do   
      
end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getOSBWE0Ind(E0Ind)
   double precision, intent(OUT)  :: E0Ind(hwlen)

   call getWinVector(Eig0, E0Ind)

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getWinVector(X, Y)
   double precision, intent(IN)   :: X(pLen(1))
   double precision, intent(OUT)  :: Y(hwLen)

   double precision, allocatable :: tX(:)
   integer :: i, ierr

   if (phwLen>0) then
       allocate(tX(phwLen))
   else
       allocate(tX(1))
   end if

   do i = 1, phwLen
      tx(i) = X(pInd(i))
   end do

   call MPI_ALLGATHERV(tX, phwLen, MPI_DOUBLE_PRECISION, Y, gCnt, sCnt,   &
                          MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)   
   deallocate(tX)

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getWinVector_CX(X, Y)
   double complex, intent(IN)   :: X(pLen(1))
   double complex, intent(OUT)  :: Y(hwLen)

   double complex, allocatable :: tX(:)
   integer :: i, ierr

   if (phwLen>0) then
       allocate(tX(phwLen))
   else
       allocate(tX(1))
   end if

   do i = 1, phwLen
      tx(i) = X(pInd(i))
   end do

   call MPI_ALLGATHERV(tX, phwLen, MPI_DOUBLE_COMPLEX, Y, gCnt, sCnt,   &
                          MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, IERR)   
   deallocate(tX)

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


