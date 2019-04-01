!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Subroutine to calculate and save E0, X0, H0 using Sinc1       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  logical function getDVR_Sinc1(isCMU, N, Xin, Mass, potVR, M, x0, E0, H0)
     implicit none
     logical, intent(IN)          :: isCMU
     integer, intent(IN)          :: N, M
     double precision, intent(IN) :: Xin,Mass
     external :: potVR
     double precision, intent(OUT):: x0(M), E0(M), H0(M, M)

     double precision,allocatable :: XR(:),HR(:,:) 
     logical :: dvr2
     integer :: info, maxNDVR

     getDVR_Sinc1=.false.;    maxNDVR=2*N+1
  
     if (maxNDVR < M) return
  
     allocate(XR(maxNDVR), HR(maxNDVR,maxNDVR), stat=info)
     if (info /= 0) return

     call Sinc1_DVR(isCMU, N, xin, mass, XR, HR)
     getDVR_Sinc1 = DVR2(maxNDVR, XR, HR, PotVR, M, X0, E0, H0)

     deallocate(XR, HR)

  end 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  logical function getDVRSF_Sinc1(isCMU,N,Xin,Mass,potVR,M,x0,E0,H0,filename)
     implicit none
     logical, intent(IN)          :: isCMU
     integer, intent(IN)          :: N, M
     double precision, intent(IN) :: Xin,Mass
     external :: potVR
     character(len=*),intent(IN)  :: filename
     double precision, intent(OUT):: x0(M), E0(M), H0(M, M)

     double precision,allocatable :: XR(:),HR(:,:),HEIG(:) 
     logical :: dvr_step2, dvr_step3
     integer :: info, maxNDVR

     getDVRSF_Sinc1=.false.;    maxNDVR=2*N+1
  
     if (maxNDVR < M) return
  
     allocate(XR(maxNDVR), HR(maxNDVR,maxNDVR), HEIG(maxNDVR), stat=info)
     if (info /= 0) return

     call Sinc1_DVR(isCMU, N, xin, mass, XR, HR)

     getDVRSF_Sinc1 = DVR_Step2(maxNDVR,XR,HR,PotVR,HEIG)
     E0(1:M)=HEIG(1:M)

     if (getDVRSF_Sinc1) then
        open(99,file=filename,status='replace',form='unformatted',iostat=info)
        if (info==0) then
           write(99) maxNDVR, -xin, xin, mass
           write(99) XR, HEIG, HR
           close(99)
           getDVRSF_Sinc1 = DVR_Step3(MAXnDVR, M, XR, HR, E0, X0, H0) 
        else
           close(99)
        end if

     end if

     deallocate(XR, HR, HEig)

  end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  logical function getDVRRF_Sinc1(M, X0, E0, H0, filename)
     implicit none
     integer, intent(IN)           :: M
     double precision, intent(OUT) :: X0(M),E0(M),H0(M,M)
     character(len=*),intent(IN)  :: filename

     double precision,allocatable :: XR(:),ER(:),HR(:,:) 
     logical :: dvr_step3
     integer :: info, Nmax
     double precision :: xmin,xmax,mass

     getDVRRF_Sinc1=.false.
    
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

     getDVRRF_Sinc1 = DVR_Step3(Nmax,M,XR,HR,E0,X0,H0)

     deallocate(XR, ER, HR)

  end 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
