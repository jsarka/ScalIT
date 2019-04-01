!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getFullViAll(Vi, Vx)

   double precision, intent(IN)  :: Vi(myVi%mEnd(sF))
   double precision, intent(OUT) :: Vx(pLen(sF))

   call getLevelViAll_Grid(sF,Vi,Vx)

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getLevelViAll(level, Vi, Vx)
   integer, intent(IN) :: level
   double precision, intent(IN) :: Vi(myVi%mEnd(level))
   double precision, intent(OUT):: Vx(pLen(level))

   if (myNode%nodNum(level)>1) then
      call getLevelViAll_Grid(level,Vi,Vx)
   else
      call getLevelViAll_Seq(level,Vi,Vx)
   end if

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getLevelViAll_Grid(level, Vi, Vx)
   integer, intent(IN) :: level
   double precision, intent(IN)  :: Vi(myVi%mEnd(level))
   double precision, intent(OUT) :: Vx(nin(level),sN(level))

   integer :: i, j, k, ind
   double precision :: tmp
   integer(kind=MPI_OFFSET_KIND) :: pos, pos1, colPos

   do i = 1, sN(level)
      pos1 = myRES%gPos(level) + (i-1)*myconf%gDim(level) 

      do j = 1, nin(level)
         pos = pos1 + j - 1  
         tmp=1.0D0

         do k = 1, level
	    colPos = pos - ((pos-1)/myconf%gLen(k))*myconf%gLen(k) 
            ind = (colPos-1)/myconf%gDim(k)             
  	    tmp = tmp * Vi(myVi%mStart(k)+ind)
         end do

         Vx(j,i) = tmp
      end do
   end do

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getLevelViAll_Seq(level, Vi, Vx)
   integer, intent(IN) :: level
   double precision, intent(IN)  :: Vi(myVi%mEnd(level))
   double precision, intent(OUT) :: Vx(plen(level))

   integer :: i, j, ind
   double precision :: tmp
   integer(kind=MPI_OFFSET_KIND) :: pos, colPos

   do i = 1, plen(level)
      pos = myRES%gPos(level) + (i-1)
      tmp=1.0D0

      do j = 1, level
	 colPos = pos - ((pos-1)/myconf%gLen(j))*myconf%gLen(j) 
         ind = (colPos-1)/myconf%gDim(j) 
  	 tmp = tmp * Vi(myVi%mStart(j)+ind)
      end do

      Vx(i) = tmp
   end do

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        Vx = V(level-1)*V(level-2)*...*V(2)*V(1)               c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getLevelVx(level, bk, row, Vi, Vx)
   integer, intent(IN) :: level,bk,row
   double precision, intent(IN)  :: Vi(myVi%mEnd(level))
   double precision, intent(OUT) :: Vx(nin(level))

   integer :: i, j, ind
   double precision :: tmp
   integer(kind=MPI_OFFSET_KIND) :: pos1, pos, colPos


   pos1 = myRES%gPos(level)+(bk-1)*myconf%gLen(level)+     &
          (row-1)*myconf%gDim(level) 
   do i = 1, nin(level)
      pos = pos1 + (i-1)
      tmp=1.0D0

      do j = 1, level-1
	 colPos = pos - ((pos-1)/myconf%gLen(j))*myconf%gLen(j) 
         ind = (colPos-1)/myconf%gDim(j) 
  	 tmp = tmp * Vi(myVi%mStart(j)+ind)
      end do

      Vx(i) = tmp
   end do

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
