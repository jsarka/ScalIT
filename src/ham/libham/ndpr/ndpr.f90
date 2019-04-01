!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Module to deal with H=T+V(r)=T+j*(j+1)/2mr^2+V(r)      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
module ndpr
   implicit none
   include 'ndpr.data'

   contains
      include 'ndpr.io'    ! read VR, and parameters
      include 'ndpr.size'  ! get DVR sizes
      include 'ndpr.dpr'   ! direct-product DVR
      include 'ndpr.xyz'   ! non-direct-product DVR
      include 'ndpr.rtp'

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function init()
         integer :: i, info, getPjmSize
         double precision :: cLege(nGL), cChev(nGC)

         init = .false.
         if (useSP) then
             print *, ' Using Spline Function to calculate 1D potential'
         else
             print *, ' Using Fitting functions to calculate 1D potential'
         end if

         allocate(nSize(lmax+1),nBase(lmax+1), nTSize(lmax+1), stat=info)
         if (info /= 0) return

         ! Original X1 and H1
         allocate(X0(nMax), H0(nMax,nMax), stat=info)
         if (info/=0) return
         call Sinc2_XH(nMax, rMin, rMax, mass, X0, H0)
 
         ! Gauss-Legendre & Gauss-Chebeshev integrations
         pjmNum = getPjmSize(lmax)
         allocate(wLege(nGL),pjmLege(nGL,pjmNum),pjmChev(nGC,pjmNum),stat=info)
         if (info /= 0) return

         call YjNodes(nGL, cLege, wLege)
         call AllYjmPolys(nGL, lmax, cLege, pjmLege)
 
         do i = 1, NGC
            cChev(i)=dcos((i-0.5D0)*PI/NGC)
         end do
         call AllYjmPolys(nGC, lmax, cChev, pjmChev)

         init = .true.

      end function
!**************************************************

      subroutine final()
         if (allocated(spR))  deallocate(spR)
         if (allocated(spV))  deallocate(spV)
         if (allocated(spM))  deallocate(spM)

         if (allocated(wLege))    deallocate(wLege)
         if (allocated(pjmLege))  deallocate(pjmLege)
         if (allocated(pjmChev))  deallocate(pjmChev)

         if (allocated(nSize))  deallocate(nSize)
         if (allocated(nTSize)) deallocate(nTSize)
         if (allocated(nBase))  deallocate(nBase)
         if (allocated(X0))     deallocate(X0)
         if (allocated(H0))     deallocate(H0)
        
      end subroutine
   
end module
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


