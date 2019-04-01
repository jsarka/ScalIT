!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Module to calculate <Fmk|phi|Fmk>, <pjm|phi|Pjm>       c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
module integralPjm
   implicit none

   double precision,parameter :: SQRT_1_2  =0.70710678118655D0 !sqrt(1/2)
   double precision,parameter :: SQRT_2_3  =0.81649658092773D0 !sqrt(2/3)
   double precision,parameter :: SQRT_4_3  =1.15470053837925D0 !sqrt(4/3)
   double precision,parameter :: SQRT_16_15=1.03279555898864D0 !sqrt(16/15)
   double precision,parameter :: SQRT_8_9  =0.94280904158206D0 !sqrt(8/9)
   double precision,parameter :: PI=3.1415926535897932D0

   private
   integer :: nGL, nGC, jmax   ! input parameters
   integer :: pjmNum
   double precision, allocatable :: wLege(:),pjmLege(:,:),pjmChev(:,:)

   public :: initIntPjm, finalIntPjm, getIntPjm, getIntPjmElem

   contains

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function initIntPjm(j0Max, nL, nC)
         integer, intent(IN) :: j0Max, nL, nC

         integer :: i, info, getPjmSize
         double precision :: cLege(nGL), cChev(nGC)

         call finalIntPjm()

         initIntPjm = .false.
         jmax=j0Max; nGL=nL; nGC=nC;  
 
         ! Gauss-Legendre & Gauss-Chebeshev integrations
         pjmNum = getPjmSize(jmax)
         allocate(wLege(nGL),pjmLege(nGL,pjmNum),pjmChev(nGC,pjmNum),stat=info)
         if (info /= 0) return

         call YjNodes(nGL, cLege, wLege)
         call AllYjmPolys(nGL, jmax, cLege, pjmLege)
 
         do i = 1, NGC
            cChev(i)=dcos((i-0.5D0)*PI/NGC)
         end do
         call AllYjmPolys(nGC, jmax, cChev, pjmChev)

         initIntPjm = .true.

      end function
!*************************************************************************

      subroutine finalIntPjm()

         if (allocated(wLege))    deallocate(wLege)

         if (allocated(pjmLege))  deallocate(pjmLege)

         if (allocated(pjmChev))  deallocate(pjmChev)
        
      end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         Calculate <Fmk|cos(phi)|Fmk>                         c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   double precision function getIntPjmElem(j1, m1, j2, m2)
      integer, intent(IN) :: j1,m1,j2,m2
      double precision :: cTh

      integer :: m12, ind1, ind2, ind0, ind
      integer :: getPjmPos

      m12  = m1 + m2
      ind1 = getPjmPos(jmax, j1, m1)  
      ind2 = getPjmPos(jmax, j2, m2)

      if (m12/2*2 == m12) then
          ind0 = getPjmPos(jmax,2, 2)
          ind  = getPjmPos(jmax,1, 0)
          cth = SUM(wLege(1:nGL)*pjmLege(1:nGL,ind1)*pjmLege(1:nGL,ind0)       &
                 *pjmLege(1:nGL,ind2))
          cth = cth * SQRT_2_3
      else
          ind0 = getPjmPos(jmax,1, 1)
          ind  = getPjmPos(jmax, 1, 0)
          cth = SUM(pjmChev(1:nGC,ind1)*pjmChev(1:nGC,ind0)*pjmChev(1:nGC,ind) &
                *pjmChev(1:nGC,ind2))
          cth = -cth * SQRT_8_9*PI/NGC
      end if 

      getIntPjmElem = cth

   end function

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         Calculate <P(jm)|sin(a).cos(a)|P(j'm')>              c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine getIntPjm(j1, m1, j2, m2, sTh, cTh)
      integer, intent(IN) :: j1,m1,j2,m2
      double precision, intent(OUT) :: sTh, cTh

      integer :: m12, ind1, ind2, ind0, ind
      integer :: getPjmPos

      m12  = m1 + m2
      ind1 = getPjmPos(jmax, j1, m1)  
      ind2 = getPjmPos(jmax, j2, m2)

      if (m12/2*2 == m12) then
         ind0 = getPjmPos(jmax, 2, 2)
         sth = SUM(pjmChev(1:nGC,ind1)*pjmChev(1:nGC,ind0)*pjmChev(1:nGC,ind2))
         sth = sth * SQRT_16_15*PI/NGC
         ind0 = getPjmPos(jmax,1, 0)
         cth = SUM(wLege(1:nGL)*pjmLege(1:nGL,ind1)*pjmLege(1:nGL,ind0)       &
                 *pjmLege(1:nGL,ind2))
         cth = cth * SQRT_2_3
     else
         ind0 = getPjmPos(jmax,1, 1)
         ind  = getPjmPos(jmax, 1, 0)
         cth = SUM(pjmChev(1:nGC,ind1)*pjmChev(1:nGC,ind0)*pjmChev(1:nGC,ind) &
                *pjmChev(1:nGC,ind2))
         cth = -cth * SQRT_8_9*PI/NGC
         ind0 = getPjmPos(jmax,1, 1)
         sth = SUM(wLege(1:nGL)*pjmLege(1:nGL,ind1)*pjmLege(1:nGL,ind0)    &
                 *pjmLege(1:nGL,ind2))
         sth = -sth * SQRT_4_3
     end if 

  end subroutine

end module
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


