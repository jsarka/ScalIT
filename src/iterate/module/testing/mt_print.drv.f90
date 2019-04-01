!ccccccccccccccccccccccccccccccccccccccccccccccc
!c        Testing program for mosbtype         c
!ccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccc
!c   comment out MPI_COMM_SPLIT function call   c
!c    at line 49 in file mosbtype.cal.h         c
!cccccccccccccccccccccccccccccccccccccccccccccccc
program mosbType_test    
   use mosbtype
   implicit none

   TYPE(CConvSimple) :: mySConv
   TYPE(CConv)     :: myConv
   TYPE(COsbw)     :: myOsbw

   print *,'  =================================='
   print *,'       Testing MOSBType: Print      '
   print *, ' =================================='
   print *
   print *, ' Read input parameters from standard input. '

   read(*,*) mySconv
   read(*,*) myConv
   read(*,*) myosbw

   print *
   print *, ' Print Simple Conv'
   call printSConv(mysconv)

   print *
   print *, ' Print Conv'
   call printCConv(myconv)

   print *
   print *, ' Print OSBW'
   call printOSBW(myosbw)

   print *, ' =========  Finish Testing ========'
   print * 

end program
!cccccccccccccccccccccccccccccccccccccccccc




