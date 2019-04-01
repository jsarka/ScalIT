!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Subroutine to read/get eigenvectors for wavefunction  c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c              Get the PIST wavefunction                       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
program test_pist
   implicit none
   character(len=128) :: fName
   integer :: Num, dN

   double precision, allocatable :: E0(:),E1(:),E2(:),V1(:,:),V2(:,:),work(:)
   integer :: M1, M2, NP, N0, info, lwork

   read(*,*) fname
   read(*,*) Num, dN

   open(99,file=fname,status='old',form='unformatted',iostat=info)
   if (info/=0) then
       print *, ' Error in open PIST file:', fName
       STOP
   end if

   read(99) NP, N0

   lwork = 3*NP
   allocate(E0(NP), E1(NP),E2(NP),V1(NP, NP), V2(NP,NP), work(lwork), stat=info)
   if (info/=0) then
      print *, ' Error in allocating Memory!'
      close(99);  STOP
   end if 

   read(99) E0(1:NP), V1(1:NP,1:NP)
   close(99)

   M1 = Min(NP, Num)
   M2 = M1 - ABS(dN)
   if (M2<0) M2 = 2

   print *
   print *, ' Current PIST Size:', NP
   print *
   print *, ' Eigen Values:'
   print *, E0

   V2(1:NP,1:NP)=V1(1:NP,1:NP)
   call DSYEV('N','U', M1, V1, NP, E1, work, lwork, info)
   
   call DSYEV('N','U', M2, V2, NP, E2, work, lwork, info)

   print * 
   write(*, 10) M1, M2
   print *, E1(1:M2)-E2(1:M2)
   deallocate(E0, E1, E2, v1, v2, work)

   10 FORMAT( ' Eigen values Difference between N1=', I6, ', and N2=', I6, '.')

end program 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
