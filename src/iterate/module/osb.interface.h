!cccccccccccccccccccccccccccccccc
!c      Public Interface        c
!cccccccccccccccccccccccccccccccc

interface readOSB
    module procedure readOSBSTD
    module procedure readOSBFile
end interface

!cccccccccccccccccccccccccccccccccc
!c          public data           c
!cccccccccccccccccccccccccccccccccc
!public :: sF, sN, sBJ, sQMR
!public :: sConv, sOSBW
!public :: sCX, sNDVR, sST, sAP, sHOSBW
!public :: sHC, sPC, sHOSB, sVOSB, sEig, sHW, sVX, sPT
!public :: myBlk, myDim, myLen
!public :: sQMRConvType

!private :: myHOSB, myH0, myVOSB, myDep
!private :: H0,H0CX,RES,OUTH,OUTHCX,DEP,DEPCX
!private :: HW,HWCX,HWR,HWRCX,HWTMP,HWTMPCX,WX,WY,WXCX,WYCX,WIPIV,WIND

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
public :: initOSB, finalOSB, initOSBW, finalOSBW
public :: readOSB, printOSB
public :: loadH0, loadRES, loadResEig, loadAP, loadDEP, loadOutH
public :: loadHOSB,saveHOSB, loadVOSBEig, saveVOSBEig
public :: loadHW, saveHW

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
public :: BlockDiag, BlockDiagUser
public :: HX,   HX0,   EHX
public :: HX_DX,HX0_DX,EHX_DX,EHX_PDX,EHX_MDX,EHX_PPDX,EHX_PRDX,EHX_MPDX,EHX_MRDX
public :: HX_CX,HX0_CX,EHX_CX,EHX_PCX,EHX_MCX,EHX_PPCX,EHX_PRCX,EHX_MPCX,EHX_MRCX
public :: PX, PX0, EPX, EPXD1, EPXD2, EPXW
public :: getOSBWSize, getOSBWIndex, getOSBWH

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
