
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c               Calculate <Fmk|phi|Fmk>                        c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   double precision function getIntFmkElem(m1, k1, m2, k2)
      integer, intent(IN) :: m1, k1, m2, k2
      double precision,parameter :: SQRT_1_2=0.70710678118655D0 !sqrt(1/2)
      double precision,parameter :: PI=3.1415926535897932D0

      double precision :: fmk

      fmk = 0.0D0
      if (k1==k2) then
         if ((m1+m2) == 0) fmk = fmk + k1 
         if ((m1-m2) == 0) fmk = fmk + 1.0D0
         if (m1==0) fmk = SQRT_1_2*fmk
         if (m2==0) fmk = SQRT_1_2*fmk
         fmk =PI*fmk
      else
         if ((m1+m2) /= 0) fmk = fmk + 1.0D0/(m1+m2)
         if ((m1-m2) /= 0) fmk = fmk + DBLE(k1)/(m2-m1)
         fmk = -fmk
         if (m1==0) fmk = SQRT_1_2*fmk
         if (m2==0) fmk = SQRT_1_2*fmk
      end if

      getIntFmkElem = fmk

   end function
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c         Calculate <F|sin(a).cos(a)|F>                        c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine getIntFmk(m1, k1, m2, k2, sfmk, cfmk)
      integer, intent(IN) :: m1, k1, m2, k2
      double precision, intent(OUT) :: sfmk, cfmk

      double precision,parameter :: SQRT_1_2=0.70710678118655D0 !sqrt(1/2)

      cfmk = 0.0D0; sfmk = 0.0D0
      if (k1/=k2) then
         if ((m2+1) == m1) sfmk = sfmk + k1 
         if ((m1+1) == m2) sfmk = sfmk + k2
         if ((m1+m2) == 1) sfmk = sfmk - 1.0D0
         sfmk = -0.5D0*sfmk
         if (m1==0) sfmk = SQRT_1_2*sfmk
         if (m2==0) sfmk = SQRT_1_2*sfmk
      else
         if ((m2+1) == m1) cfmk = cfmk + 1.0D0
         if ((m1+1) == m2) cfmk = cfmk + 1.0D0
         if ((m1+m2) == 1) cfmk = cfmk + k1
         cfmk = 0.5D0*cfmk
         if (m1==0) cfmk = SQRT_1_2*cfmk
         if (m2==0) cfmk = SQRT_1_2*cfmk
      end if

   end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
