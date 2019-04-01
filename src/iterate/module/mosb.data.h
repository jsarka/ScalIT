!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Define the data used for MOSB package           c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   integer, parameter :: STDFH = 5
   integer, parameter :: MAX_FNAME = 128

           !sJOB, job type: bound state, resonance state, CRP
   integer, parameter :: JOB_BOUND=1,JOB_RES1=2, JOB_RES2=3,   &
                         JOB_CRP=4,  JOB_CRP1=5, JOB_CRP2=6

           ! various preconditioners: sHC, sPC, defalut: (H-E0)*X
           ! sHC, sPC: H*X,(H-E0)*X,(H+/-iAP)*X
   integer, parameter :: TA0=0, TA1=1, TAAP0=2, TMAP0=3

           ! sHC, sPC: (H-E0+iAP[APR,APP])*X
   integer, parameter :: TAAP=100, TAAPP=101, TAAPR=102

           !sHC, sPC: (H-E0-iAP[APR,APP])*X
   integer, parameter :: TMAP=200, TMAPP=201, TMAPR=202

           !sOSB, type of preconditioner: OSB, OSBD and OSBW
   integer, parameter :: TOSB=1,TOSBD1=2,TOSBD2=3,TOSBW=4

   integer, parameter :: NTYPE_X=1,NTYPE_DX=2,NTYPE_CX=3

   integer, parameter :: rootID = 0

   logical, parameter :: myDebug=.false.

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   integer :: sF, sN(FMAX)   

   TYPE(CConvSimple) :: sBJ, sQMR
   TYPE(CConv)       :: sConv
   TYPE(COsbw)       :: sOSBW

   integer :: sJOB, sOSB
   logical :: sCX,  sNDVR, sST, sAP 
   integer :: sHC,  sPC,  sHij
   integer :: sHOSB, sVOSB, sHW, sVX,  sPT
   logical :: sDEP(FMAX) 

   character(len=MAX_FNAME) :: fH0, fRES
   character(len=MAX_FNAME) :: fAPP,  fAPR, fOUTH, fDEP
   character(len=MAX_FNAME) :: fHOSB, fVOSB, fEig
   character(len=MAX_FNAME) :: fHW, fVX, fPT
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   integer :: id, nNodes
   integer :: dbSize, cxSize
   logical :: totalDep, srMode=.true.
   integer :: nin(FMAX), nout(FMAX), blk(FMAX), plen(FMAX)
   integer :: sNMax, pmax
    
   TYPE(GDataInfo) :: myconf      ! information for global data
   TYPE(MNodeInfo) :: myNode      ! node information
   TYPE(MDataInfo) :: myData      ! information for local data 
   TYPE(CSeqInfo)  :: myH0, myVi
   TYPE(MGridInfo) :: myHOSB, myVOSB, myRES, myDep

   double precision, allocatable :: RES(:), EIG0(:), ResSeq(:), VOSB(:)
   double precision, allocatable :: APP(:), APR(:),  AP(:)
   double precision, allocatable :: H0(:),  HOSB(:),  OUTH(:),  DEP(:)
   double complex, allocatable   :: H0CX(:),HOSBCX(:),OUTHCX(:),DEPCX(:)

   integer :: sQMRConvType=-1    ! conv. and initial choice for QMR
   integer :: sConvType=1        ! conv. testing for PIST/Lan
   integer :: sCXConvType=0      ! conv. testing for complex value.
   integer :: sPISTConvType=0
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




