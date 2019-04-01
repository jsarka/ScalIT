!
! Testing for CG and get3index 
!
PROGRAM TEST_cg
   implicit none

   integer :: jmax(3)
   INTEGER, allocatable :: jbase(:),jlen0(:), jLen1(:),jLen2(:),jLen3(:), jLen4(:)
   integer, allocatable :: jInd0(:,:),jInd1(:,:), jInd2(:,:),jInd3(:,:),jInd4(:,:)
   integer :: getCgJMSize, getCgCoeffSize, getCgCoeffSize1
   integer :: jmSize, jfSize, jffSize
   logical :: s1
   integer ::  I, OPT = 1

   DO WHILE (OPT /= 0)
       PRINT *, 'Input jmax for Cg symbol: j1, j2, j3 '
       read  *, jmax(1:3)
       jmSize = getCgJMSize(jmax) 
       jfsize = getCgCoeffSize(jmax)
       print *, 'Size of non-zero Cg coefficients. Total:',jmSize, jfSize

       print *, 'Whether to Show Length and base: T=Show, F=Not Show'
       read *, s1
       if (s1) then
          allocate(jLen0(jmSize), jBase(jmSize))
          jffSize = getCgCoeffSize1(jmax,jmSize, jLen0, jBase)
          print *
          print *, 'Size:', jffSize, sum(jLen0)
          print *, 'JLength:', jLen0
          print *, 'JBase:', jBase
          deallocate(jLen0, jBase)
       end if

       print *
       print *, 'Whether to Show index: T=Show, F=Not Show'
       read *, s1
       if (s1) then
           allocate(jLen1(jmSize), jLen2(jmSize),                       &
                    jLen3(jmSize), jLen4(jmsize), jInd0(4,jmSize),      &
                    jInd1(4,jmSize), jInd2(4,jmSize), jInd3(4, jmSize), &
                    jInd4(4,jmSize))
           print *

           print *, ' index<->(jm) for 3j symbols:'
           print *, 'Index   J   M   New Pos    Error'
             
           call getCgLen   (jmax, jmSize, JLen1)
           call getCgLen0  (jmax, jmSize, JLen2)
          print *, 'jLen2', jLen2         
           call getCgIndex (jmax, jmSize, jInd1)
           call getCgIndex0(jmax, jmSize, jind2)
           call getCgLenInd(jmax, jmSize, jLen3, jInd3)
           call getCgLenInd0(jmax, jmSize, jLen4,jind4)

           print *, 'Length:'
           print *, 'Label, Length, Index:'
           do i = 1, jmSize
              print *
              print *, i, jLen1(i), jInd1(:,i)
              print *, i, jLen3(i), jInd3(:,i)
              print *, i, jLen2(i), jind2(:,i)
              print *, i, jLen4(i), jInd4(:,i)
           end do

           deallocate(jLen1,jLen2,jLen3,jLen4,jInd0,jInd1, jInd2, jInd3, jInd4)
       end if
       print *
       PRINT *, 'INPUT OPTION for loop: O FOR EXIT'
       READ *, OPT
    END DO
END
