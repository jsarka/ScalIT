!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  Calculate Hij: use more memory, but less computing time      c
!c    All return .false. if memory allocation occurs error       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function getHOSBH0(M, ind1, ind2, nH)
   integer, intent(IN) :: M,ind2(M)
   integer(kind=MPI_OFFSET_KIND),intent(IN) :: ind1(M)
   double precision, intent(OUT) :: nH(M,M)

   double precision :: nH0(M,M)
   double precision, allocatable :: Vij(:,:)
   double precision, allocatable :: allVi(:),allVj(:)

   integer :: i, j, k, level, mylevel
   integer :: info, ierr

   getHOSBH0 = .false.
   allocate(Vij(myVi%mLen,M), stat=info)
   if (info/=0) return

   do i=1, M
      call getFullVi( i, Vij(1,i) )
   end do

   nH0(1:M,1:M)=0.0D0
   do i = 1, M
      if (ind2(i)/=0) nH0(i,i) = Eig0(ind2(i))
   end do

   do i = 1, M
      do j = i+1, M         
  	 do k = sF, 1, -1
	    if ((blkInd(k,i)==blkInd(k,j)).AND.(sNInd(k,i)/=sNInd(k,j))) then
	         level = k; exit
            end if
         end do
 
         if (grpInd(level,i)==myNode%grpID(level))   then

            allocate(allVi(nin(level)), allVj(nin(level)),stat=info)
            if (info/=0) return

            call getLevelVx(level,blkInd(level,i),sNInd(level,i),Vij(1,i),allVi)

            call getLevelVx(level,blkInd(level,j),sNInd(level,j),Vij(1,j),allVj)

            nH0(i,j) = HOSB_dotProd(blk(level),nin(level),sN(level),      &
                       blkInd(level,i),sNInd(level,i),sNInd(level,j),     &
                       allVi, allVj, HOSB(myHOSB%pStart(level)))  

            nH0(j,i) = nH0(i,j) 

            deallocate(allVi, allVj) 

         end if      
      end do    ! end of j
   end do       ! end of i

   call MPI_AllReduce(nH0, nH, M*M, MPI_DOUBLE_PRECISION, MPI_SUM,  &
               MPI_COMM_WORLD, ierr)

   deallocate(Vij)

   getHOSBH0 = .true.

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function getFOSBH0(M, ind1, ind2, nH)
   integer, intent(IN) :: M,ind2(M)
   integer(kind=MPI_OFFSET_KIND),intent(IN) :: ind1(M)
   double precision, intent(OUT) :: nH(M,M)

   double precision :: nH0(M,M)
   double precision, allocatable :: Vij(:,:)
   double precision, allocatable :: allVi(:),allVj(:)

   integer :: i, j, k, mylevel, level
   integer :: fh, ierr, info

   getFOSBH0 = .false.
   allocate(Vij(myVi%mLen,M), stat=info)
   if (info/=0) return

   do i=1, M
      call getFullVi( i, Vij(1,i) )
   end do

   nH0(1:M,1:M)=0.0D0
   do i = 1, M
      if (ind2(i)/=0) nH0(i,i) = Eig0(ind2(i))
   end do

   call MPI_File_Open(MPI_COMM_WORLD,fHOSB,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)

   do level = sF, 1, -1

      allocate(HOSB(myHOSB%pSize(level)),allVi(nin(level)),allVj(nin(level)),&
               stat=info)
      if (info/=0) return

      if (myNode%nodNum(level)>1) then
         call MReadDataGrid(fh, dbSize,myHOSB%gPos(level),myconf%gDim(level),&
               myData%pDim(level),sN(level),HOSB(myHOSB%pStart(level)),ierr)
      else
         call MReadData(fh, dbSize,myHOSB%gPos(level),myHOSB%pSize(level),   &
                   HOSB(myHOSB%pStart(level)), ierr)
      end if

      do i = 1, M

         do j = i+1, M         

  	    do k = sF, 1, -1
	       if ((blkInd(k,i)==blkInd(k,j)).AND.(sNInd(k,i)/=sNInd(k,j))) then
	          myLevel = k; exit
               end if
            end do
 
	    if ( level == myLevel ) then

               if (grpInd(level,i)==myNode%grpID(level))   then

             	   call getLevelVx(level,blkInd(level,i),sNInd(level,i), &
                                   Vij(1,i),allVi)

                   call getLevelVx(level,blkInd(level,j),sNInd(level,j), &
                                   Vij(1,j),allVj)

                   nH0(i,j) = HOSB_dotProd(blk(level),nin(level),sN(level), &
                            blkInd(level,i),sNInd(level,i),sNInd(level,j),  &
                            allVi, allVj, HOSB(myHOSB%pStart(level)))

                   nH0(j,i) = nH0(i,j) 

               end if
            end if        
         end do
      end do               

      deallocate(HOSB, allVi, allVj)

   end do

   call MPI_FILE_CLOSE(fh, ierr)

   call MPI_AllReduce(nH0, nH, M*M, MPI_DOUBLE_PRECISION, MPI_SUM,  &
               MPI_COMM_WORLD, ierr)

   deallocate(Vij)

   getFOSBH0 = .true.

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
logical function getVOSBH0(M, ind1, ind2, nH)
   integer, intent(IN) :: M,ind2(M)
   integer(kind=MPI_OFFSET_KIND),intent(IN) :: ind1(M)
   double precision, intent(OUT) :: nH(M, M)

   integer :: i, j, ierr, info
   double precision :: nH0(M, M)


   double precision, allocatable :: Vij(:,:)
   double precision, allocatable :: allVi(:,:),allVj(:,:),HVj(:,:)

   getVOSBH0 = .false.
   allocate(Vij(myVi%mLen,M), allVi(myData%pDim(sF),sN(sF)),  &
            allVj(myData%pDim(sF),sN(sF)),HVj(myData%pDim(sF),&
            sN(sF)), stat=info)
   if (info/=0) return

   do i=1, M
      call getFullVi( i, Vij(1,i) )
   end do

   nH0(1:M,1:M)=0.0D0
   do i = 1, M
      if (ind2(i)/=0) nH0(i,i) = Eig0(ind2(i))
   end do

   do i = 1, M

      call getFullViAll(Vij(1,i),allVi)

      do j = i+1, M

          call getFullViAll(Vij(1,j),allVj)

          call HijX(myData%pDim(sF)*sN(sF),allVj, HVj)

          nH0(i,j) = MA_dotProd(myData%pDim(sF)*sN(sF), allVi, HVj)

          nH0(j,i) = nH0(i,j)

      end do      

   end do

   call MPI_AllReduce(nH0, nH, M*M, MPI_DOUBLE_PRECISION, MPI_SUM, &
                      MPI_COMM_WORLD, ierr)

   deallocate(Vij, allVi, allVj, HVj)

   getVOSBH0 = .true.

end function
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


