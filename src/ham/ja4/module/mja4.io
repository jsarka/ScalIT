!*************************************************
!*   Read Input parameters and distribute them   *
!*************************************************
!****************************************************
subroutine DistributeInput()
   integer, parameter  :: nInt = 5
   integer, parameter  :: nLog = 3

   logical :: logData(nLog)
   integer :: intData(nInt), info
   
   if (myid==rootID) then
      logData(1)=parity; logData(2)=useSP
   end if
   call MPI_BCAST(logData, nLog, MPI_LOGICAL,rootID,MPI_COMM_WORLD, ierr)
   parity = logData(1);  useSP=logData(2)
   
   if (myid==rootID) then
      intData(1)=JTol; 
      intData(2)=FcFlag;  intData(3)=CbFlag
      intData(4)=ReFlag;  intData(5)=absFlag;
   end if
   call MPI_BCAST(intData,nInt, MPI_INTEGER, rootID, MPI_COMM_WORLD, ierr)
   JTol=intData(1) ; 
   FcFlag=intData(2);  CbFlag=intData(3)
   ReFlag=intData(4);  absFlag=intData(5)

   call MPI_BCAST(jmax, NA, MPI_INTEGER, rootID, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(NDVR, NR+1, MPI_INTEGER, rootID, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(NGI,  NA, MPI_INTEGER, rootID, MPI_COMM_WORLD, ierr)

   call MPI_BCAST(MASS,NR,MPI_DOUBLE_PRECISION,rootID,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(RE,  NR,MPI_DOUBLE_PRECISION,rootID,MPI_COMM_WORLD,ierr)
 
   call MPI_BCAST(fH0,FILENAMELEN,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(fH0GM,FILENAMELEN,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(fRe,FILENAMELEN,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)

   if (useSP) then
      call MPI_BCAST(fSPVRlr1,FILENAMELEN,MPI_CHARACTER,rootID,MPI_COMM_WORLD, &
                     ierr)
      call MPI_BCAST(fSPVRlr2,FILENAMELEN,MPI_CHARACTER,rootID,MPI_COMM_WORLD, &
                     ierr)
      call MPI_BCAST(fSPVRBR, FILENAMELEN,MPI_CHARACTER,rootID,MPI_COMM_WORLD, &
                     ierr)
   end if

   if ((absFlag==ABS_ONE) .OR. (absFlag==ABS_TWO))   &
      call MPI_BCAST(fABS,FILENAMELEN,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)

end subroutine
!*******************************************************************************


!**********************************************************************
subroutine DistributeDVRData()

   call MPI_BCAST(spNum,NR,MPI_INTEGER,rootID,MPI_COMM_WORLD,ierr)

   if (myid/=rootID) then
      allocate(splr1(spnum(1)),spvlr1(spNum(1)),spMlr1(spNum(1)), &
               splr2(spnum(2)),spvlr2(spNum(2)),spMlr2(spNum(2)), &
               spBR(spnum(3)),spvBR(spNum(3)),spMBR(spNum(3)))
   end if

   Call MPI_BCAST(splr1, spNum(1), MPI_DOUBLE_PRECISION, rootID,  &
                    MPI_COMM_WORLD, ierr)
   Call MPI_BCAST(spVlr1,spNum(1), MPI_DOUBLE_PRECISION, rootID,  &
                    MPI_COMM_WORLD, ierr)
   Call MPI_BCAST(splr2, spNum(2), MPI_DOUBLE_PRECISION, rootID,  &
                    MPI_COMM_WORLD, ierr)
   Call MPI_BCAST(spVlr2,spNum(2), MPI_DOUBLE_PRECISION, rootID,  &
                    MPI_COMM_WORLD, ierr)
   Call MPI_BCAST(spBR, spNum(3), MPI_DOUBLE_PRECISION, rootID,  &
                    MPI_COMM_WORLD, ierr)
   Call MPI_BCAST(spVBR,spNum(3), MPI_DOUBLE_PRECISION, rootID,  &
                    MPI_COMM_WORLD, ierr)

   if (myid/=rootID) then
      call spline(spNum(1),splr1,spVlr1,spyp1(1),spyp2(1),spMlr1)
      call spline(spNum(2),splr2,spVlr2,spyp1(2),spyp2(2),spMlr2)
      call spline(spNum(3),spBR,spVBR,spyp1(3),spyp2(3),spMBR)
   end if

end subroutine
!**********************************************************************

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine distributeR_VR()

    Call MPI_BCAST(lr1,   NDVR(1), MPI_DOUBLE_PRECISION, workID(1),  &
                    MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(vlr1,  NDVR(1), MPI_DOUBLE_PRECISION, workID(1),  &
                    MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(momlr1,NDVR(1), MPI_DOUBLE_PRECISION, workID(1),  &
                    MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(Elr1,  NDVR(1), MPI_DOUBLE_PRECISION, workID(1),  &
                    MPI_COMM_WORLD, ierr)

    Call MPI_BCAST(lr2,   NDVR(2), MPI_DOUBLE_PRECISION, workID(1),  &
                    MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(vlr2,  NDVR(2), MPI_DOUBLE_PRECISION, workID(1),  &
                    MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(momlr2,NDVR(2), MPI_DOUBLE_PRECISION, workID(1),  &
                    MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(Elr2,  NDVR(2), MPI_DOUBLE_PRECISION, workID(1),  &
                    MPI_COMM_WORLD, ierr)

    Call MPI_BCAST(Br,    NDVR(3), MPI_DOUBLE_PRECISION, workID(1),  &
                    MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(vBR,   NDVR(3), MPI_DOUBLE_PRECISION, workID(1),  &
                    MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(momBR, NDVR(3), MPI_DOUBLE_PRECISION, workID(1),  &
                    MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(EBR, NDVR(3), MPI_DOUBLE_PRECISION, workID(1),    &
                    MPI_COMM_WORLD, ierr)
end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

