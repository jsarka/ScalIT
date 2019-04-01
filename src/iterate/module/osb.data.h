!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Define the data used for OSB package           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c               Constants for various work                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           !sOSB, type of preconditioner: OSB, OSBD and OSBW
   integer, parameter :: TOSB=1,TOSBD1=2,TOSBD2=3,TOSBW=4   

           !sJOB, job type: bound state, resonance state, CRP
   integer, parameter :: JOB_BOUND=1,JOB_RES1=2, JOB_RES2=3,  &
                         JOB_CRP=4,  JOB_CRP1=5, JOB_CRP2=6

           !various preconditioners: sHC, sPC, defalut: (H-E0)*X
           !sHC, sPC: H*X,(H-E0)*X,(H+/-iAP)*X
   integer, parameter :: TA0=0, TA1=1 ,TAAP0=2, TMAP0=3  
           !sHC, sPC: (H-E0+iAP[APR,APP])*X
   integer, parameter :: TAAP=100,TAAPP=101,TAAPR=102 
           !sHC, sPC: (H-E0-iAP[APP,APR])*X
   integer, parameter :: TMAP=200,TMAPP=201,TMAPR=202 

           ! [-SEQDIR, SEQDIR]=direct access
   integer, parameter :: SEQDIRNUM=100  

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   integer :: sF, sN(FMAX) 
   logical :: sDEP(FMAX)  
   integer :: sOSB,sJOB
   logical :: sCX, sNDVR, sST, sAP     
   TYPE(CConvSimple) :: sBJ, sQMR
   TYPE(CConv)       :: sConv
   TYPE(COsbw)       :: sOSBW
   integer :: sHOSB, sVOSB, sHW, sVX, sPT   ! sHOSB=1:100 direct save,>100:seq
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   character(len=MAX_FNAME) :: H0File, RESFile
   character(len=MAX_FNAME) :: APPFile, APRFile, OUTHFile, DEPFile
   character(len=MAX_FNAME) :: HOSBFile, VOSBFile, EigFile
   character(len=MAX_FNAME) :: HWFile, VXFile, PTFile
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


   integer :: sQMRConvType=-1    ! conv. and initial choice for QMR
   integer :: sConvType=1        ! conv. testing for PIST/Lan
   integer :: sCXConvType=0      ! conv. testing for complex value.
   integer :: sPISTConvType=0

   integer :: myBlk(FMAX), myDim(FMAX), myClen(FMAX)
   integer :: myLen, outLen
   integer (kind=8) :: hwLen    
   logical :: totalDep

   integer :: sHC, sPC, sHij     ! for H*X, P*X and Hij

   logical :: srMode=.TRUE.

   TYPE(CDataInfo) :: myVi, myH0, myHOSB, myVOSB, myDep

   double precision, allocatable :: RES(:),EIG0(:),APP(:),APR(:),SQ_APR(:),AP(:),VOSB(:)
   double precision, allocatable :: H0(:), HOSB(:),OUTH(:),DEP(:)
   double complex, allocatable   :: H0CX(:),HOSBCX(:),OUTHCX(:),DEPCX(:)

   double precision, allocatable :: HW(:,:),  HWR(:,:),  HWTMP(:,:)
   double precision, allocatable :: WX(:),WY(:)
   double complex, allocatable   :: HWCX(:,:),HWRCX(:,:),HWTMPCX(:,:)
   double complex, allocatable   :: WXCX(:),WYCX(:)
   integer, allocatable  :: WIPIV(:), WIND(:)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

