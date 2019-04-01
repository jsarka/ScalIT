!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine DistributeInput()
   integer, parameter  :: nInt = 7 
   integer, parameter  :: nLog = 3
! The following line was changed from nDbl value 5 to 7 to avoid out of
!   bounds warnings when compiling.  Up to now, we didn't use Absorption
!   potentials which seems to be broken here.  
!   - Corey Petty 3/2/2015
   integer, parameter  :: nDbl = 7

   logical :: logData(nLog)
   integer :: intData(nInt), info
   double precision :: dbData(nDbl)

   if (myid==rootID) then
      logData(1)=parity;  logData(2)=useSP;  logData(3)=saveMode
   end if
   call MPI_BCAST(logData,nLog,MPI_LOGICAL,rootID,MPI_COMM_WORLD,ierr)
   parity = logData(1);  useSP=logData(2);   !saveMode=.TRUE.

   if (myid==rootID) then
      intData(1)=JTol;   intData(2)=FcFlag;  intData(3)=CbFlag
      intData(4)=ReFlag; intData(5)=absFlag; 
      intData(6)=en(1);  intData(7)=en(2)
   end if
   call MPI_BCAST(intData,nInt,MPI_INTEGER,rootID,MPI_COMM_WORLD,ierr)
   JTol=intData(1);   FcFlag=intData(2);  CbFlag=intData(3)
   ReFlag=intData(4); absFlag=intData(5);
   en(1)=intData(6);  en(2)=intData(7)

   if (myid==rootID) then
      dbData(1)=Ecutoff; 
      dbData(2)=A0(1);    dbData(3)=A0(2)
      dbData(4)=Rabs0(1); dbData(5)=Rabs0(2); 
      dbData(6)=Rabs1(1); dbData(7)=Rabs1(2)
   end if
   call MPI_BCAST(dbData,nDbl,MPI_DOUBLE_PRECISION,rootID,MPI_COMM_WORLD,ierr)
   Ecutoff=dbData(1);
   A0(1)=dbData(2);     A0(2)=dbData(3)
   Rabs0(1)=dbData(4);  Rabs0(2)=dbData(5);
   Rabs1(1)=dbData(6);  Rabs1(2)=dbData(7)

   call MPI_BCAST(jmax, NA, MPI_INTEGER, rootID, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(NDVR, NR+1, MPI_INTEGER, rootID, MPI_COMM_WORLD, ierr)
   call MPI_BCAST(NGI,  NA, MPI_INTEGER, rootID, MPI_COMM_WORLD, ierr)

   call MPI_BCAST(MASS,NR,MPI_DOUBLE_PRECISION,rootID,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(RE, NR,MPI_DOUBLE_PRECISION,rootID,MPI_COMM_WORLD,ierr)

   call MPI_BCAST(fH0,FILENAMELEN,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(fH0GM,FILENAMELEN,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)
   if (ReFlag/=0)  &
      call MPI_BCAST(fRe,FILENAMELEN,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(fVRlr,FILENAMELEN,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(fVRBR,FILENAMELEN,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)

   if (useSP) then
      call MPI_BCAST(fSpVRlr,FILENAMELEN,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(fSpVRBR,FILENAMELEN,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)
   end if

   if ((absFlag==ABS_ONE) .OR. (absFlag==ABS_TWO))   &
      call MPI_BCAST(fABS,FILENAMELEN,MPI_CHARACTER,rootID,MPI_COMM_WORLD,ierr)

end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine DistributeDVRData()

   call MPI_BCAST(spNum,NR,MPI_INTEGER,rootID,MPI_COMM_WORLD,ierr)
   if (myid/=rootID) then
      allocate(splr(spnum(1)),spvlr(spNum(1)),spMlr(spNum(1)), &
               spBR(spnum(2)),spvBR(spNum(2)),spMBR(spNum(2)) )
   end if

   Call MPI_BCAST(splr, spNum(1), MPI_DOUBLE_PRECISION, rootID,  &
                  MPI_COMM_WORLD, ierr)
   Call MPI_BCAST(spVlr,spNum(1), MPI_DOUBLE_PRECISION, rootID,  &
                  MPI_COMM_WORLD, ierr)

   Call MPI_BCAST(spBR, spNum(2), MPI_DOUBLE_PRECISION, rootID,  &
                  MPI_COMM_WORLD, ierr)
   Call MPI_BCAST(spVBR,spNum(2), MPI_DOUBLE_PRECISION, rootID,  &
                  MPI_COMM_WORLD, ierr)

   if (myid/=rootID) then
      call spline(spNum(1),splr,spVlr,spyp1(1),spyp2(1),spMlr)
      call spline(spNum(2),spBR,spVBR,spyp1(2),spyp2(2),spMBR)
   end if
          
end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine distributeR_VR()

    Call MPI_BCAST(lr,  NDVR(1), MPI_DOUBLE_PRECISION, workID(1),  &
                    MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(vlr, NDVR(1), MPI_DOUBLE_PRECISION, workID(1),  &
                    MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(momlr,NDVR(1), MPI_DOUBLE_PRECISION, workID(1), &
                    MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(Elr, NDVR(1), MPI_DOUBLE_PRECISION, workID(1),  &
                    MPI_COMM_WORLD, ierr)

    Call MPI_BCAST(Br,  NDVR(2), MPI_DOUBLE_PRECISION, workID(1),  &
                    MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(vBR, NDVR(2), MPI_DOUBLE_PRECISION, workID(1),  &
                    MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(momBR,NDVR(2), MPI_DOUBLE_PRECISION, workID(1), &
                    MPI_COMM_WORLD, ierr)
    Call MPI_BCAST(EBR, NDVR(2), MPI_DOUBLE_PRECISION, workID(1),  &
                    MPI_COMM_WORLD, ierr)
end subroutine
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
