!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!c           Hij part of OSB:                                     c
!c   Calculate OSBW matrix elements for the states |Ei-E0|<DE     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                Calculate Hij using HOSB                        c
!c        Here all HOSB should be stored in the memory            c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getOSBWHCX(nSize, ind, nH)
   integer (kind=8), intent(IN)          :: nSize
   integer, intent(IN)          :: ind(nSize)
   double complex, intent(OUT) :: nH(nSize, nSize)

   if (sST) then
       if (.NOT. getHOSBH0_CX(nSize, ind, nH))   &
          print *, ' Error in calculating Window H0.'
       deallocate(HOSBCX)
   else
      if (sHOSB/=0) then
         if (.NOT.getFOSBH0_CX(nSize, ind, nH))  &
            call getVOSBH0_CX(nSize, ind, nH)  
      else
          call getVOSBH0_CX(nSize, ind, nH)      
      end if
   end if
     
end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  logical function getHOSBH0_CX(cnt, ind, Hij)
      integer (kind=8), intent(IN)  :: cnt
      integer, intent(IN)           :: ind(cnt)
      double complex, intent(OUT) :: Hij(cnt, cnt) 

      integer :: I, J, K, level, info
      integer :: bkInd(sF,cnt), colInd(sF,cnt) 

      getHOSBH0_CX = .false.
      if ((.NOT.sST) .OR. (.NOT. allocated(HOSB)) .OR.   &
          (size(HOSB) /= myHOSB%mLen) )   return

      call getLevelIndex(cnt, ind, bkInd, colInd)
 
      do I = 1, cnt
         do J = I + 1, cnt
 	    do K = sF, 1, -1
	       if ( (bkInd(k,i)==bkInd(k,j)) .AND.       &
                    (colInd(k,i)/=colInd(k,j)) ) then
	          level = K;  exit
               end if
            end do

            Hij(i,j) = getHij_CX(level,bkInd(1:sF,i), colInd(1:sF,i), &
                    colInd(1:sF,j), HOSBCX(myHOSB%mStart(level)))

            Hij(j, i) = CONJG(Hij(i, j))
         end do
         Hij(i, i) =  EIG0(ind(i)) 
      end do

      getHOSBH0_CX = .true.

 end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  logical function getFOSBH0_CX(cnt, ind, Hij)
      integer (kind=8), intent(IN)  :: cnt
      integer, intent(IN)           :: ind(cnt)
      double complex, intent(OUT) :: Hij(cnt, cnt) 

      integer :: I, J, k, level, info, num, sumOff
      integer :: bkInd(sF,cnt), colInd(sF,cnt)
      double complex :: dbcx
 
      getFOSBH0_CX = .false.
      call getLevelIndex(cnt, ind, bkInd, colInd)

      inquire(IOLENGTH=num) dbcx
      open(99,FILE=HOSBFILE,STATUS='OLD',FORM='UNFORMATTED',       &
              access='direct',recl=num,IOSTAT=info)
      if (info /=0) return 

      do level = 1, sF
         allocate(HOSBCX(myHOSB%mSize(level)),stat=info)
         if (info/=0) then
             close(99); return
         end if

         ! read(99,IOSTAT=info) HOSBCX(1:myHOSB%mSize(level))
         ! if (info/=0) then
         !    close(99); return
         ! end if

         sumOff = (Sum(sN(1:level))-sN(level))*myLen
  	 do i = 1, sN(level)*myLen
 	    num = sumOff + i
 	    write(99,REC=num) HOSBCX(i)
         end do

         do I = 1, cnt
            do J = I + 1, cnt
	       if ( (bkInd(level,i)==bkInd(level,j)) .AND.          &
                    (colInd(level,i)/=colInd(level,j)) )   then
                    Hij(i,j) = getHij_CX(level,bkInd(1:sF,i),colInd(1:sF,i),  &
                           colInd(1:sF,j), HOSBCX(myHOSB%mStart(level)))
                    Hij(j, i) = CONJG(Hij(i, j))
               end if
            end do
         end do
         deallocate(HOSBCX)
      end do
      close(99)

      do i = 1, cnt
         Hij(i, i) =  EIG0(ind(i)) 
      end do
     
      getFOSBH0_CX = .true.

 end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        get the whole Hij for states within the window          c
!c                     only use the VOSB                          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getVOSBH0_CX(cnt, ind, Hij)
   integer (kind=8), intent(IN)  :: cnt
   integer, intent(IN)           :: ind(cnt)
   double complex, intent(OUT) :: Hij(cnt, cnt)

   integer :: I, J
   double precision :: X0(myLen), X1(myLen)
   double complex   :: Y0(myLen), Y1(mylen), tmp(mylen)
   integer :: bkInd(sF,cnt), colInd(sF,cnt)

   call getLevelIndex(cnt, ind, bkInd, colInd)

   do I = 1, cnt
      call getVi(bkInd(1:sF,i),colInd(1:sF,i),X0)   ! cpu_time=T
      Y0(1:mylen) = X0(1:mylen)
      do J = I + 1, cnt
         call getVi(bkInd(1:sF,j),colInd(1:sF,j),X1) 

         Y1(1:myLen) = X1(1:myLen)
         call HijX_CX(mylen, Y1, tmp)               ! cpu_time=10T

         Hij(I, J) = dot_product(Y0(1:myLen), tmp(1:mylen))
         Hij(J,I)  = CONJG(Hij(I,J))     

      end do

      Hij(I, I) =  EIG0(ind(I)) 
   end do

end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Calculate the specified HIJ = V^T|H|V  via HOSB        c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
double complex function getHij_CX(level, bk, row, col, HB)
      integer, intent(IN) :: level, bk(sF), row(sF), col(sF)
      double complex, intent(IN) :: HB(myDim(level),sN(level),   &
                                         sN(level),myBlk(level))

      double precision :: VI(myDim(level)), VJ(myDim(level))
      double complex   :: Y0(myDim(level))

      call getLevelVi(level-1, bk, row, myDim(level), vi) 
      call getLevelVi(level-1, bk, col, myDim(level), vj)

      Y0(1:myDim(level)) = Vi(1:myDim(level))*Vj(1:myDim(level))
      getHij_CX=sum(Y0(1:myDim(level))*HB(1:myDim(level),        &
                    row(level),col(level),bk(level)))

              !DOTVTDV(len, vi, HB, vj)     
end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
