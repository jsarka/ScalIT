!
! Module to store Gaunt coefficients and 
! extract them from the memory
!
module threejmod
  implicit none
  private    
  
  public :: init3j, final3j, get3jCoeff0, get3jCoeff1, get3jCoeff2
  public :: getCGCoeff0, getCGCoeff1,getCGCoeff2

  integer :: jMax3j(3) = (/0,0,0/)
  integer :: jmNum3j = 0, cfNum3j=0, jmmax3j=0
  integer, allocatable :: jmLen3j(:), jmBase3j(:), sqrtj3(:)
  double precision, allocatable :: coeff3j(:)
  logical :: st = .FALSE.


  contains

!*******************************************************************
  logical function init3j(j0)
     integer, intent(IN) :: j0(3)  
     
     integer :: info, i
     integer :: get3jJmSize, get3jCoeffSize1
   
     init3j = .FALSE.;       call final3j()
     jMax3j(1:3) = j0(1:3)

     jmNum3j = get3jJmSize(jmax3j)
     jmmax3j = max(jmax3j(1),jmax3j(2),jmax3j(3)) + 1
     allocate(jmLen3j(jmNum3j), jmBase3j(jmNum3j), sqrtj3(jmmax3j), stat=info)
     if (info /= 0)  return
     do i = 0, jmmax3j
        sqrtj3(i) = DSQRT(1.0D0+i+i)
     end do

     cfNum3j = get3jCoeffSize1(jmax3j, jmNum3j, jmLen3j, jmBase3j)
     allocate(coeff3j(cfNum3j), stat=info)
     if (info /= 0)  return

     call fill3j(jmax3j, jmNum3j, jmLen3j, jmBase3j, cfNum3j, coeff3j)
     init3j = .TRUE.
     st   = .TRUE.

  end function

!******************************************************************
  subroutine final3j()
      
      if (allocated(jmLen3j))   deallocate(jmLen3j)
      if (allocated(jmBase3j))  deallocate(jmBase3j)
      if (allocated(sqrtj3))    deallocate(sqrtj3)
      if (allocated(coeff3j))   deallocate(coeff3j)

      st = .FALSE.

  end subroutine

!*****************************************************************
  double precision function get3jCoeff0(jm)
       integer, intent(IN) :: jm(3,2)
       
       integer :: pos, get3jPos
       logical :: isNonZero3j
      
       get3jCoeff0 = 0.0D0
       if (ST)  then
          if ( isNonZero3j(jm(1,1),jm(1,2),jm(2,1),jm(2,2),jm(3,1),jm(3,2))) then
             pos = get3jPos(jmNum3j, jmLen3j, jmBase3j, jm)
             if ((pos == 0) .OR. (pos > cfNum3j)) then
                print *, ' Error to get Clebsh Gordon Coefficients at ',  &
                          jm, ' with position:', pos
             else                
                if (pos<0) then
                   get3jCoeff0 = -coeff3j(-pos)
                else
                   get3jCoeff0 = coeff3j(pos)
                end if
             end if
          end if
       else
          print *, ' Please call init() before using this function'
       end if

  end function

  double precision function get3jCoeff1(j, m)
       integer, intent(IN) :: j(3), m(3)

       integer :: jm(3,2)

       jm(1:3,1) = j(1:3); jm(1:3,2)=m(1:3)
       get3jCoeff1 = get3jCoeff0(jm)
  end function

  double precision function get3jCoeff2(j1, m1, j2, m2, j3, m3)
       integer, intent(IN) :: j1, m1, j2, m2, j3, m3

       integer :: jm(3,2)
       
       jm(1,1)=j1; jm(2,1)=j2; jm(3,1)=j3
       jm(1,2)=m1; jm(2,2)=m2; jm(3,2)=m3
       get3jCoeff2 = get3jCoeff0(jm)

  end function

!*****************************************************************
  double precision function getCGCoeff0(jm)
       integer, intent(IN) :: jm(3,2)
       
       integer :: pos, getCGPos, jmi
       logical :: isNonZeroCG
      
       getCGCoeff0 = 0.0D0
       if (ST)  then
          if ( isNonZeroCG(jm(1,1),jm(1,2),jm(2,1),jm(2,2),jm(3,1),jm(3,2))) then
             pos = getCGPos(jmNum3j, jmLen3j, jmBase3j, jm)
             if ((pos ==0) .OR. (pos > cfNum3j)) then
                print *, ' Error to get Clebsh Gordon Coefficients at ',  &
                          jm, ' with position:', pos
             else
                if (pos < 0) then 
                   getCGCoeff0 = -coeff3j(-pos)
                else
                   getCGCoeff0 = coeff3j(pos)
                end if
                getCGCoeff0 = SQRTJ3(jm(3,1)+1)*getCGCoeff0
                jmi = jm(1,1)-jm(2,1)+jm(3,2)
                if (jmi/2*2 /= jmi) getCGCoeff0 = -getCGCoeff0
             end if
          end if
       else
          print *, ' Please call init() before using this function'
       end if

  end function

  double precision function getCGCoeff1(j, m)
       integer, intent(IN) :: j(3), m(3)

       integer :: jm(3,2)

       jm(1:3,1) = j(1:3); jm(1:3,2)=m(1:3)
       getCGCoeff1 = getCGCoeff0(jm)
  end function

  double precision function getCGCoeff2(j1, m1, j2, m2, j3, m3)
       integer, intent(IN) :: j1, m1, j2, m2, j3, m3

       integer :: jm(3,2)
       
       jm(1,1)=j1; jm(2,1)=j2; jm(3,1)=j3
       jm(1,2)=m1; jm(2,2)=m2; jm(3,2)=m3
       getCGCoeff2 = getCGCoeff0(jm)

  end function

!*****************************************************************


end module


