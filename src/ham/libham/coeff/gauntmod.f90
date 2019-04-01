!
! Module to store Gaunt coefficients and 
! extract them from the memory
!
module gauntmod
  implicit none
  private
 
  public :: initGA, finalGA, getGauntCoeff0, getGauntCoeff1,getGauntCoeff2

  integer :: jMaxga(3) = (/0,0,0/)
  integer :: jmNumga = 0, cfNumga=0
  integer, allocatable :: jmLenga(:), jmBasega(:)
  double precision, allocatable :: coeffga(:)
  logical :: stga = .FALSE.

  contains

!*******************************************************************
logical function initGA(j0)
     integer, intent(IN) :: j0(3)  
     
     integer :: info
     integer :: getGauntJmSize, getGauntCoeffSize1
   
     initGA = .FALSE.;       call finalGA()
     jMaxga(1:3) = j0(1:3)

     jmNumga = getGauntJmSize(jmaxga)
     allocate(jmLenga(jmNumga), jmBasega(jmNumga), stat=info) 
     if (info /= 0)  return

     cfNumga = getGauntCoeffSize1(jmaxga, jmNumga, jmLenga, jmBasega)
     allocate(coeffga(cfNumga), stat=info)
     if (info /= 0)  return

     call fillGaunt(jmaxga,jmNumga,jmLenga,jmBasega,cfNumga,coeffga)
     initGA = .TRUE.
     stga   = .TRUE.
  end function

!******************************************************************
  subroutine finalGA()

      if (allocated(jmLenga))  deallocate(jmLenga)
      if (allocated(jmBasega)) deallocate(jmBasega)
      if (allocated(coeffga))  deallocate(coeffga)
      stga = .FALSE.

  end subroutine

!*****************************************************************
  double precision function getGauntCoeff0(jm)
       integer, intent(IN) :: jm(3,2)
       
       integer :: pos, getGauntPos
       logical :: isNonZeroGaunt
      
       getGauntCoeff0 = 0.0D0
       if (STga)  then
          if ( isNonZeroGaunt(jm(1,1),jm(1,2),jm(2,1),jm(2,2),    &
                              jm(3,1),jm(3,2))) then
             pos = getGauntPos(jmNumga, jmLenga, jmBasega, jm)
             if ((pos < 1) .OR. (pos > cfNumga)) then
                print *, ' Error to get Gaunt Coefficients at ',  &
                          jm, ' with position:', pos
             else
                getGauntCoeff0 = coeffga(pos)
             end if
          end if
       else
          print *, ' Please call init() before using this function'
       end if

  end function

  double precision function getGauntCoeff1(j, m)
       integer, intent(IN) :: j(3), m(3)

       integer :: jm(3,2)

       jm(1:3,1) = j(1:3); jm(1:3,2)=m(1:3)
       getGauntCoeff1 = getGauntCoeff0(jm)
  end function

  double precision function getGauntCoeff2(j1, m1, j2, m2, j3, m3)
       integer, intent(IN) :: j1, m1, j2, m2, j3, m3

       integer :: jm(3,2)
       
       jm(1,1)=j1; jm(2,1)=j2; jm(3,1)=j3
       jm(1,2)=m1; jm(2,2)=m2; jm(3,2)=m3
       getGauntCoeff2 = getGauntCoeff0(jm)

  end function

!*****************************************************************
 

end module


