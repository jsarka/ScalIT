!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getMOSBWH(M, ind1, ind2, nH)
     integer, intent(IN) :: M,ind2(M)
   integer(kind=MPI_OFFSET_KIND),intent(IN) :: ind1(M)
   double precision, intent(OUT) :: nH(M,M)

   if ( sST ) then 
      if (id==rootID)  print *, ' Call getHOSBH0 to calculate Hij of OSBW'
      if (.NOT.getHOSBH0(M, ind1, ind2, nH))      &
          call getHOSBH1(M, ind1, ind2, nH)
   else    
      if (sHOSB<0) then
         if (id==rootID)  print *, ' Call getFOSBH0 to calculate Hij of OSBW'

         if (.NOT.getFOSBH0(M, ind1, ind2, nH))   &
             call getFOSBH1(M, ind1, ind2, nH) 
      else
         if (id==rootID)  print *, ' Call getVOSBH0 to calculate Hij of OSBW'
         if (.NOT.getVOSBH0(M, ind1, ind2, nH))   &
             call getVOSBH1(M, ind1, ind2, nH)  
      end if
   end if

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getMOSBWHCX(M, ind1, ind2, nH)
   integer, intent(IN) :: M, ind2(M)
   integer(kind=MPI_OFFSET_KIND),intent(IN) :: ind1(M)
   double complex, intent(OUT) :: nH(M, M)

   if (sST) then
      if (id==rootID)  print *, ' Call getHOSBH0 to calculate Hij of OSBW'
      if (.NOT. getHOSBH0CX(M, ind1, ind2, nH))      &
         call getHOSBH1CX(M, ind1, ind2, nH)
   else
      if (sHOSB<0) then
         if (id==rootID)  print *, ' Call getFOSBH0 to calculate Hij of OSBW'
         if (.NOT. getFOSBH0CX(M, ind1, ind2, nH))   &
             call getFOSBH1CX(M, ind1, ind2, nH)
      else
         if (id==rootID)  print *, ' Call getVOSBH0 to calculate Hij of OSBW'
         if (.NOT. getVOSBH0CX(M, ind1, ind2, nH))   &
           call getVOSBH1CX(M, ind1, ind2, nH)  
      end if
   end if

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

