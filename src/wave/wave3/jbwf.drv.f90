!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Subroutines to calculate wave-function        c
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
program test_jbwf
   implicit none
   character(len=*),parameter :: fVlr='~/work2/ne30/input/ne3lr.dat'
   character(len=*),parameter :: fVBr='~/work2/ne30/input/ne3Br.dat'
   character(len=*),parameter :: fVTH='~/work2/ne30/input/RE.dat'
   character(len=*),parameter :: fPT='~/work2/ne30/input/pt.dat'
   
   integer, parameter :: JTol=0, JMAX=60, NSTATE=4
   logical, parameter :: parity=.true.
 
   integer :: NMAX(3), NS(3), ntotal, p1NMax,NSIND(NSTATE)

   double precision, allocatable :: pjm(:), vpjm(:), Rlr(:), RBR(:), RTH(:)
   double precision, allocatable :: VBR(:,:),Vlr(:,:),VTh(:,:),VP(:,:),  &
                       V0(:,:,:),V1(:)   
   integer, allocatable :: jkInd(:), jIndex(:), kIndex(:)

   integer :: getPjmPos, getPjmSize, get3Size
   logical :: getPartDVRV_Sinc2,getPISTV_Index, getConV

   double precision :: BRMin, BRMax, lrMin, lrMax, dBR, dLR
   double precision :: BR(3), lr(3), cth(3)
   integer :: i, j, k, i0, j0, k0, info
   double precision :: P0(3,NSTATE)

   print *, ' Input JB parameters: BR, lr, cos(theta)'
   read(*,*) BR(1),lr(1),cth(1)
   
   NMAX(1)=6000; NMAX(2)=6000
   NS(1)=40;  NS(2)=60
   p1NMax = getPjmSize(JMAX); NMAX(3)=get3Size(parity,JTol,JMAX)
   NS(3) = NMAX(3); ntotal=NS(1)*NS(2)*NS(3)

   allocate( pjm(p1NMax), vpjm(NMax(3)), Rlr(NS(1)), RBR(NS(2)), RTH(NS(3)), &
             jkInd(NMax(3)), jIndex(NMax(3)), kIndex(NMax(3)))

   if (NMAX(3)==NS(3)) then
      allocate(Vlr(NMax(1),NS(1)), VBR(NMax(2),NS(2)), VP(NTotal, NState), &
               V0(NS(1),NS(2),NS(3)), V1(Ntotal),stat=info)
   else
      allocate(Vlr(NMax(1),NS(1)), VBR(NMax(2),NS(2)), VP(NTotal, NState), &
               VTh(NMax(3),NS(3)), V0(NS(1),NS(2),NS(3)), V1(Ntotal),stat=info)
   end if
   if (info/=0) stop
 
   call get3Index(parity, JTol, Jmax, NMax(3), jIndex, kIndex)
   do i = 1, Nmax(3)
      jkInd(i)=getPjmPos(Jmax,jIndex(i),kIndex(i))
   end do

   do i = 1, NSTATE
      NSInd(i)=i
   end do

  if ((.NOT.(getPartDVRV_Sinc2(fVlr,NMax(1),NS(1),LRMin,LRMax,Vlr))) .OR. &
       (.NOT.(getPartDVRV_Sinc2(fVBR,NMax(2),NS(2),BRMin,BRMax,VBR))) .OR. &
       (.NOT.(getPISTV_Index(fPT,Ntotal,NState,NSInd,VP))))  then
       if (allocated(VTh)) deallocate(VTh)
       deallocate(Vlr,VBR,VP,V0, V1)
       deallocate( pjm, vpjm, Rlr, RBR, RTH, jkInd, jIndex, kIndex)
       stop
   end if

   if (NMAX(3)/=NS(3)) then
      if (.NOT.getConV(fVTh, NMax(3), NS(3), VTH)) then
         deallocate(Vlr,VBR,VTh,VP,V0,V1); 
         deallocate( pjm, vpjm, Rlr, RBR, RTH, jkInd, jIndex, kIndex)
         stop
      end if
   end if
            
   dLR = (lrMax-LRMin)/NMax(1);    dBR = (BRMax-BRMin)/NMax(2)
   call JBECH(BR(1),lr(1),cth(1),BR(2),lr(2),cth(2),BR(3),lr(3),cth(3))

   print *, ' JB parameters:'
   print *, '  BR:        lr:        cos(theta)'
   do i = 1, 3
      print *, BR(i), lr(i), cth(i)
   end do

   do i = 1, 3
      call AllYjmPolys(1,jmax, cTh(i), pjm)

      vpjm(1:NMax(3))=pjm(jkInd(1:NMax(3)))

      if (Nmax(3)==NS(3)) then
         RTh(1:NS(3))=vpjm(1:NMax(3))
      else
       ! fpjm(1,NS(3))=vpjm(1,Ntheta)*VTH(Ntheta,NS(3))
        do k = 1, NS(3)
           RTH(k) = dot_product(vpjm(1:NMax(3)),VTh(1:NMax(3),k))
        end do
      end if

      k = (BR(i)-BRMin)/dBR+1
      if (k<1) k=1
      if (k>NMax(2)) k=NMax(2)
      RBR(1:NS(2)) = VBR(k,1:NS(2))
 
      k = (lr(i)-lrMin)/dlr+1
      if (k<1) k=1
      if (k>NMax(1)) k=NMax(1)
      Rlr(1:NS(1)) = Vlr(k,1:NS(1))

      do k0 = 1, NS(3)
         do j0 = 1, NS(2)
            do i0 = 1, NS(1)
               v0(i0, j0, k0) = Rlr(i0) * RBr(j0) * Rth(k0)
            end do
         end do 
      end do

      V1=reshape(V0,(/Ntotal/))
      do k =1, Nstate
         P0(k,I)=dot_product(VP(1:Ntotal,k),V1(1:Ntotal)) 
      end do

   end do

   print *
   print *, ' Wave function for equival  values:'
   do i = 1, NSTATE
      print *, i, P0(i,1), P0(i,2), P0(i,3)
   end do

   if (allocated(VTh)) deallocate(VTh)
   deallocate(Vlr,VBR,VP,V0, V1)
   deallocate( pjm, vpjm, Rlr, RBR, RTH, jkInd, jIndex, kIndex)

   
end program

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
