!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Subroutine to calculate and save E0, X0, H0 using Sinc2       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  logical function getDVR_Sinc2(isCMU, N, A0, B0, Mass, potVR, M, x0, E0, H0)
     implicit none
     logical, intent(IN)          :: isCMU
     integer, intent(IN)          :: N, M
     double precision, intent(IN) :: A0, B0, Mass
     external :: potVR
     double precision, intent(OUT):: X0(M), E0(M), H0(M, M)

     double precision,allocatable :: XR(:),HR(:,:)
     logical :: dvr2
     integer :: info

     getDVR_Sinc2=.false.
     if ( N < M) return

     allocate(XR(N), HR(N,N), stat=info)
     if (info /= 0) return

     call Sinc2_DVR(isCMU, N, A0, B0, mass, XR, HR)

     getDVR_Sinc2 = DVR2(N, XR, HR, PotVR, M, x0, E0, H0) 

     deallocate(XR, HR)

  end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  logical function getDVRSF_Sinc2(isCMU,N,A0,B0,Mass,potVR,M,x0,E0,H0,filename)
     implicit none
     logical, intent(IN)          :: isCMU
     integer, intent(IN)          :: N, M
     double precision, intent(IN) :: A0, B0, Mass
     external :: potVR
     character(len=*),intent(IN)  :: filename
     double precision, intent(OUT):: X0(M), E0(M), H0(M, M)

     double precision,allocatable :: XR(:),HR(:,:), HEIG(:)
     logical :: dvr_step2, dvr_step3
     integer :: info

     getDVRSF_Sinc2=.false.
     if ( N < M) return

     allocate(XR(N), HR(N,N), HEIG(N), stat=info)
     if (info /= 0) return

     call Sinc2_DVR(isCMU, N, A0, B0, mass, XR, HR)

     getDVRSF_Sinc2 = DVR_Step2(N,XR,HR,PotVR,HEIG)
     E0(1:M)=HEIG(1:M)

     if (getDVRSF_Sinc2) then
        open(99,file=filename,status='replace',form='unformatted',iostat=info)
        if (info==0) then
           write(99) N, min(A0,B0), max(A0,B0), mass
           write(99) XR, HEIG, HR
           close(99)
           getDVRSF_Sinc2 = DVR_Step3(N, M, XR, HR, E0, X0, H0)
        else
           close(99)
        end if

     end if
     
     deallocate(XR, HR, HEIG)

  end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  logical function getDVRRF_Sinc2(M, X0, E0, H0, filename)
     implicit none
     integer, intent(IN)           :: M
     double precision, intent(OUT) :: X0(M),E0(M),H0(M,M)
     character(len=*),intent(IN)  :: filename

     double precision,allocatable :: XR(:),ER(:),HR(:,:)
     logical :: dvr_step3
     integer :: info, Nmax
     double precision :: xmin,xmax,mass

     getDVRRF_Sinc2=.false.

!     print *, 'DVR file:', filename

     open(99,file=filename,status='old',form='unformatted',iostat=info)
     if (info/=0) return

     read(99) Nmax, xmin, xmax, mass

     if (Nmax<M) then
        close(99); return
     end if

     allocate(XR(Nmax), ER(Nmax), HR(Nmax,M), stat=info)
     if (info /= 0) then
        close(99); return
     end if

     read(99) XR(1:Nmax),ER(1:Nmax),HR(1:Nmax,1:M)
     close(99)

     E0(1:M) = ER(1:M)
     getDVRRF_Sinc2 = DVR_Step3(Nmax,M,XR,HR,E0,X0,H0)

     deallocate(XR, ER, HR)

  end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!*********************************************************************

