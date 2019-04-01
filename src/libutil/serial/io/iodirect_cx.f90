!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       General subroutines to perform I/O operations:              c
!c            save data, read data, and print data                   c
!c                  Using direct access method                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                     Save data                               c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  logical function saveDataDir_CX(N, data1, fileName)
      implicit none
      integer, intent(IN) :: N
      double complex, intent(IN) :: data1(N)
      character(LEN=*), intent(IN) ::  fileName

      integer :: info
      integer(kind=8) :: num

      saveDataDir_CX = .false.

      inquire(IOLENGTH=num) data1    

      open(99,FILE=fileName, STATUS='REPLACE', FORM='UNFORMATTED', &
               ACCESS='direct',RECL=num, IOSTAT=info)
      if (info /= 0) return

      write(99, rec=1) data1

      close(99)

      saveDataDir_CX = .true. 

  end 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                     Save data2                              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  logical function saveData2Dir_CX(N, data1, data2, fileName)
      implicit none
      integer, intent(IN) :: N
      double complex, intent(IN) :: data1(N), data2(N)
      character(LEN=*), intent(IN) :: fileName

      integer :: info
      integer(kind=8) :: num
 
      saveData2Dir_CX = .false.

      inquire(IOLENGTH=num) data1  

      open(99,FILE=fileName, STATUS='REPLACE', FORM='UNFORMATTED', &
              ACCESS='direct', RECL=num, IOSTAT=info)        
      if (info/=0) return

      write(99,rec=1) Data1
      write(99,rec=2) Data2

      close(99)

      saveData2Dir_CX = .true.

   end
 
!*****************************************************************
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                     load data                               c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  logical function loadDataDir_CX(N, data1, fileName)
      implicit none
      integer, intent(IN) :: N
      double complex, intent(OUT) :: data1(N)
      character(LEN=*), intent(IN) :: fileName

      integer :: info
      integer(kind=8) :: num

      loadDataDir_CX = .false.

      inquire(IOLENGTH=num) data1  

      open(99,FILE=fileName, STATUS='OLD', FORM='UNFORMATTED', &
            ACCESS='direct', RECL=num, IOSTAT=info)
      if (info /= 0) return

      read(99, rec=1) Data1

      close(99)

      loadDataDir_CX = .true. 
  end 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                     Save data2                              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  logical function loadData2Dir_CX(N, data1, data2, fileName)
      implicit none
      integer, intent(IN) :: N
      double complex, intent(OUT) :: data1(N), data2(N)
      character(LEN=*), intent(IN) ::  fileName

      integer :: info
      integer(kind=8) :: num
 
      loadData2Dir_CX = .false.

      inquire(IOLENGTH=num) data1 

      open(99,FILE=fileName, STATUS='OLD', FORM='UNFORMATTED', &
             ACCESS='direct',IOSTAT=info)

      read(99, rec=1) Data1
      read(99, rec=2) Data2

      close(99)

      loadData2Dir_CX = .true.

   end
 
!*****************************************************************
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
