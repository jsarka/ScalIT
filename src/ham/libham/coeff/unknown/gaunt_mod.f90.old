!
! Module to store Gaunt coefficients and 
! extract them from the memory
!
module GauntMod
  implicit none
  private    ! public
  integer :: jMax(3) = (/0,0,0/)
  integer :: jmNum = 0, cfNum=0
  integer, allocatable :: jmLen(:), jmBase(:)
  double precision, allocatable :: coeff(:)
  logical :: st = .FALSE.

!  interface getGaunt
!     Module Procedure getGauntCoeff0
!     Module Procedure getGauntCoeff1
!     Module Procedure getGauntCoeff2
!  end interface

!  interface getJmNum
!     Module Procedure getJmNum0
!     Module Procedure getJMNum1
!  end interface
!
!  interface getCFNum
!     Module Procedure getCFNum0
!     Module Procedure getCFNum1
!  end interface

  contains

!*******************************************************************
  logical function initGA(j0)
     integer, intent(IN) :: j0(3)  
     
     integer :: info
     integer :: getGauntJmSize, getGauntCoeffSize
   
     initGA = .FALSE.;       call finalGA()
     jMax(1:3) = j0(1:3)

     jmNum = getGauntJmSize(jmax)
     allocate(jmLen(jmNum), jmBase(jmNum), stat=info) 
     if (info /= 0)  return

     cfNum = getGauntCoeffSize(jmax, jmNum, jmLen, jmBase)
     allocate(coeff(cfNum), stat=info)
     if (info /= 0)  return

     call fillGaunt(jmax, jmNum, jmLen, jmBase, cfNum, coeff)
     initGA = .TRUE.
     st   = .TRUE.
  end function

!******************************************************************
  subroutine finalGA()

      if (allocated(jmLen))  deallocate(jmLen)
      if (allocated(jmBase)) deallocate(jmBase)
      if (allocated(coeff))  deallocate(coeff)
      st = .FALSE.

  end subroutine

!*****************************************************************
  double precision function getGauntCoeff0(jm)
       integer, intent(IN) :: jm(3,2)
       
       integer :: pos, getGauntPos
       logical :: isNonZeroGaunt
      
       getGauntCoeff0 = 0.0D0
       if (ST)  then
          if ( isNonZeroGaunt(jm)) then
             pos = getGauntPos(jmNum, jmLen, jmBase, jm)
             if ((pos < 1) .OR. (pos > cfNum)) then
                print *, ' Error to get Gaunt Coefficients at ',  &
                          jm, ' with position:', pos
             else
                getGauntCoeff0 = coeff(pos)
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
  integer function getJ1Max()
      getJ1Max = Jmax(1)
  end function 

  integer function getJ2Max()
      getJ2Max = Jmax(2)
  end function

  integer function getJ3Max()
      getJ3Max = Jmax(3)
  end function

  integer function getJmNum0()
      getJmNum0 = jmNum
  end function

  integer function getJmNum1(jm)
      integer, intent(IN) :: jm(3)

      integer :: getGauntJMSize
     
      getJmNum1 = getGauntJMSize(jm)
  end function

  integer function getCFNum0()
      getCFNum0 = cfNum
  end function

  integer function getCFNum1(jm)
      integer, intent(IN) :: jm(3)

      integer :: info, jSize, getGauntJMSize
      integer, allocatable :: jLen(:)

      getCFNum1 = -1
      jSize = getGauntJMSize(jm)
      allocate (jLen(jSize), stat=info)
      if (info /= 0)  return
      call getGauntLen(jm, jSize, jLen)  
      getCFNum1 = sum(jLen(1:jSize))
      deallocate(jLen)
  end function

  logical function getStatus()
      getStatus = ST
  end function 

end module


