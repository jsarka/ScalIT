!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Using Sinc1(2) DVR Function as the basis function to create a   c
!c general basis function                                          c
!c                                                                 c
!c Using spline function or external function of the potential and c
!c kinetic energy from Sinc1(2) calculates and stores eigenvalues  c
!c (eigenvectors) of the Hamiltonian matrix.                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
module PreSinc
    implicit none
    private    
!c
    integer, parameter :: FN_LEN=128  !..Length of file name
    double precision, parameter :: spyp1=3.0D36, spyp2=3.0D36 !..Parameters for spline 
!c
    integer :: NMAX    !..Maximum number of SINC DVR functions
    integer :: N       !..Number of computed VBR energy eigenstates       
    integer :: spNum   !..Number of spline data
    integer :: nType   !..1: sinc1 DVR, 2:sinc2 DVR
    logical :: useSP   ! .TRUE., use Spline function, otherwise use fitting function
    double precision :: MASS         !..mass
    double precision :: xmin, xmax  
!c..x range, [-xmax,xmax] for Sinc1, [xmin,xmax] for Sinc2
!c                     
    character(len=FN_LEN) :: spFile, outFile
!c..Spline data are always in ASCII format: It has the following format:
!c..nSP, (spR(i), spVR(i), i=1, nSP)
!c
    double precision, allocatable :: spR(:), spV(:), spM(:)    
    double precision, allocatable :: H0(:),  X0(:), Eig0(:)
    integer :: MaxNDVR    !..Actual number of SINC DVR function
!c
    public :: run
!c
  contains
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine run(fitV)    !..User only needs to call this function
        external :: fitV
        double precision :: ct1, ct2, ct3, ct4
        integer :: info
        logical :: suc
!c       
        print *
        print *, '========================================'
        print *, '     Prepare Data for Matrix-Elements   '
        print *, '========================================'
!c        
        print *, ' Read Input Parameters from Standard IO'
!c
        call CPU_Time(ct1)
        call readPS()   !..Read input parameters
        call printPS()  !..Print related parameters
!c
        if (useSP) then
          if ( .NOT. readSplineData()) then             
              print *, ' Error in Read spline Data from :', spFile            
              return
           end if
        end if
!c        
        allocate(X0(maxNDVR),H0(maxNDVR*maxNDVR),Eig0(maxNDVR),stat=info)
        if (info/=0) then        
           print *, ' Error in allocating memory for H0, X0'
           return
        else
           print *, ' Perform PSO-DVR ...... '
        end if
!c
        if (useSP) then
            suc = getData(calV)
        else
            suc = getData(fitV)
        end if
!c
        print *
        if (.NOT. suc) then          
           print *, ' Error in performing PSO-DVR for H0, X0'
           return
        else       
           print *, ' Eigenvalues of the interested states:'
           print 999, Eig0(1:N)
        end if
!c
        call CPU_Time(ct3)
        call saveData()
        call CPU_Time(ct4)
        print *
        print *, ' CPU Time to save step 1st data files:', ct4-ct3
        print *
        call deallocMemory()
        call CPU_Time(ct2)
        print *
        print *, '   CPU Time for this preparation:', ct2-ct1
        print *        
        print *, '====      Finish the calculation         ==='
        print *, '============================================='
!c        
 999 format (E25.15,2x,E25.15,2x,E25.15)
    end subroutine
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine readPS()   !..Read input parameters
!c
        read(*,*) nType, mass, Nmax, N, useSP
!c
        if (nType==1) then
            read(*,*) xmax
            xmin=-xmax
            maxNDVR = 2*NMax+1
        else
            read(*,*) Xmin, xmax
            maxNDVR = NMax
        end if
!c
        read(*,'(A)') outFile
        if (useSP) read(*,'(A)') spFile
!c
    end subroutine
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
    logical function readSplineData()
       integer :: i, info
!c       
       readSplineData = .false.
       open(99,FILE=spFile, status='old', iostat=info)
!c
       if (info/=0) return
       read(99,*) spNum
!c
       allocate(spR(spNum),spV(spNum), spM(spNum), stat=info)
       if (info/=0) return
!c       
       do i=1,spNum
          read(99,*) spR(i), spV(i)
       end do
       close(99)
!c
       call spline(spNum,spR,spV,spyp1,spyp2,spM)
!c.....Given arrays spR and spV of order spNum and given values spyp1 and 
!c.....spyp2 for the first derivative of the interpolating function at points
!c.....1 and spNum, respectively, this routine returns an array spM of length spNum
!c.....which contains the second derivatives of the interpolating function.

       readSplineData = .TRUE.
!c
    end function
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine printPS()
!c
        print *
        print *, '======================================='
        if (nType==1) then
           print *, ' Using Sinc1 DVR as Basis Functions'
           print *, ' X range:[',-Xmax, ',', xmax,']'
        else
           print *, ' Using Sinc2 DVR as Basis Functions'
           print *, ' X range:[',Xmin, ',', xmax,']'
        end if
!c
        print *
        print *, ' Max. Numer of original Function:', maxNDVR
        print *, ' Number of interested States:', N
        print *, ' Mass:', mass
!c
        print *
        if (useSP) then
           print *, ' Using Splining Function to calculate 1D potential'
           print *, ' Read spline data from file: ', spFile
        else
           print *, ' Using Fitting Function to calculate 1D potential'
        end if
        print *, ' Final Results are stored in file:', outFile
!c
    end subroutine
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine calV(N, R, VR)
       integer, intent(IN) :: N
       double precision, intent(IN)  :: R(N)
       double precision, intent(OUT) :: VR(N)
!c       
       call splint(spNum, spR, spV, spM, N, R, VR)
!c..Using spNum, spR,spV and spline output spM, this routine returns a 
!c..cubic-spline interpolated value VR.
!c
    end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
    logical function getData(fitV)
       external :: fitV       
       logical :: DVR_Step2
!c
       if (nType==1) then
           call Sinc1_DVR(.false.,NMAX, xmax,mass,X0,H0) !..General SINC (-infinity:infinity)
       else
           call Sinc2_DVR(.false.,NMAX,xmin,xmax,mass,X0,H0) !..Radial SINC (0:infinity)
       end if
       getData = DVR_Step2(MaxNDVR,X0,H0,fitV,Eig0)
!c..DVR_Step2 returns eigenvlaues Eig0/ eigenvectors of (H0+Eig0).
!c
    end function
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine saveData()
!c 
        open(99,file=outFile,status='replace',form='unformatted')
        write(99) maxNDVR,xmin,xmax,mass 
        write(99) X0, Eig0, H0      
        close(99)
!c
    end subroutine
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
    subroutine deallocMemory()
!c
       if (allocated(X0))   deallocate(X0)
       if (allocated(H0))   deallocate(H0)
       if (allocated(Eig0)) deallocate(Eig0)
       if (allocated(spR))  deallocate(spR)
       if (allocated(spV))  deallocate(spV)
       if (allocated(spM))  deallocate(spM)
!c
    end subroutine
!c    
end module
!c
