!ccccccccccccccccccccccccccccccccccccccccccccccc
!c        Testing program for mosbtype         c
!ccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccc
!c   comment out MPI_COMM_SPLIT function call   c
!c    at line 49 in file mosbtype.cal.h         c
!cccccccccccccccccccccccccccccccccccccccccccccccc
program mosb_seq
   use mosbtype
   implicit none

   integer :: opt, i, j, k

   logical :: sNDVR
   integer :: nNodes, sF, sN(FMAX)
   Type(GDataInfo) :: myGd
   TYPE(CSeqInfo)  :: myH0, myVi   

   print *,'  =================================='
   print *,'        Testing MOSBType            '
   print *, ' =================================='
   print *
   print *, ' Read input parameters from standard input. '

   read(*,*) nNodes, sF
   read(*,*) sN(1:sF)
  
   call calGData(sF,sN,myGD)

   print *
   print *, ' Global Information'
   call printGdata(myGD)

   print *
   print *, '-------------------------------'
   print *
   print *, ' Vi indices'
   call calViSeq(sF,sN,myVi)
   call printSeqInfo(myVi)

   print *
   print *, '-------------------------------'
   print *
   sNDVR = .true.
   print *, ' H0 indices, sNDVR=',sNDVR
   call calH0Seq(sNDVR, sF,sN,myH0) 
   call printSeqInfo(myH0)

   print *
   print *, '-------------------------------'
   print *
   sNDVR = .false.
   print *, ' H0 indices, sNDVR=',sNDVR

   call calH0Seq(sNDVR, sF,sN,myH0) 
   call printSeqInfo(myH0)
   print *, '*********   Finish Testing  *******'

end program
!cccccccccccccccccccccccccccccccccccccccccc



