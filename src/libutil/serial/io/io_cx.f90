!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       General subroutines to perform I/O operations:              c
!c            save data, read data, and print data                   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                     Save data                               c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  logical function saveData_CX(N, data1, binMode, fileName)
      implicit none
      integer, intent(IN) :: N
      double precision, intent(IN) :: data1(N)
      logical, intent(IN) :: binMode
      character(LEN=*), intent(IN) ::  fileName

      integer :: info

      saveData_CX = .false.

      if (binMode) then
          open(99,FILE=fileName, STATUS='REPLACE', FORM='UNFORMATTED', &
               IOSTAT=info)
          if (info /= 0) return
          write(99) Data1  !(1:N)
      else
          open(99, FILE=FILENAME,STATUS='REPLACE', IOSTAT=info)
          if (info /= 0) return
          write(99, *) DATA1  !(1:N)
      end if

      close(99)
      saveData_CX = .true. 
  end 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                     Save data2                              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  logical function saveData2_CX(N, data1, data2, binMode, fileName)
      implicit none
      integer, intent(IN) :: N
      double precision, intent(IN) :: data1(N), data2(N)
      logical, intent(IN) :: binMode
      character(LEN=*), intent(IN) ::  fileName

      integer :: info
 
      saveData2_CX = .false.
      if (binMode) then
          open(99,FILE=fileName, STATUS='REPLACE', FORM='UNFORMATTED', &
               IOSTAT=info)
          if (info/=0) return
          write(99) Data1 ! (1:N)
          write(99) Data2 ! (1:N)
      else
          open(99, FILE=FILENAME,STATUS='REPLACE', IOSTAT=info)
          if (info/=0) return
          write(99, *) data1  ! (1:N)
          write(99, *) data2  ! (1:N)
      end if

      close(99)
      saveData2_CX = .true.
   end
 
!*****************************************************************
   logical function save2Data_CX(N, data1, data2, binMode, fileName)
      implicit none
      integer, intent(IN) :: N
      double complex,intent(IN) :: data1(N), data2(N)
      logical, intent(IN) :: binMode
      character(LEN=*), intent(IN) ::  fileName

      integer :: i, info

      save2Data_CX = .false.
      if (binMode) then
          open(99,FILE=fileName, STATUS='REPLACE', FORM='UNFORMATTED',  &
               IOSTAT=info)
          if (info/=0) return
          do I = 1, N
             write(99) Data1(I), data2(I)
          end do
      else
          open(99, FILE=FILENAME,STATUS='REPLACE', IOSTAT=info)
          if (info/=0) return
          do I = 1, N
             write(99, *) data1(I), data2(I)
          end do
      end if
         
      close(99)
      save2Data_CX = .true.
    end 

!*************************************************************


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                     load data                               c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  logical function loadData_CX(N, data1, binMode, fileName)
      implicit none
      integer, intent(IN) :: N
      double complex, intent(OUT) :: data1(N)
      logical, intent(IN) :: binMode
      character(LEN=*), intent(IN) ::  fileName

      integer :: info

      loadData_CX = .false.

      if (binMode) then
          open(99,FILE=fileName, STATUS='OLD', FORM='UNFORMATTED', &
               IOSTAT=info)
          if (info /= 0) return
          read(99) Data1  ! (1:N)
      else
          open(99, FILE=FILENAME,STATUS='OLD', IOSTAT=info)
          if (info /= 0) return
          read(99, *) DATA1   ! (1:N)
      end if

      close(99)
      loadData_CX = .true. 
  end 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                     load data2                              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  logical function loadData2_CX(N, data1, data2, binMode, fileName)
      implicit none
      integer, intent(IN) :: N
      double complex, intent(OUT) :: data1(N), data2(N)
      logical, intent(IN) :: binMode
      character(LEN=*), intent(IN) ::  fileName

      integer :: info
 
      loadData2_CX = .false.
      if (binMode) then
          open(99,FILE=fileName, STATUS='OLD', FORM='UNFORMATTED', &
               IOSTAT=info)
          if (info/=0) return
          read(99) Data1  !(1:N)
          read(99) Data2  !(1:N)
      else
          open(99, FILE=FILENAME,STATUS='OLD', IOSTAT=info)
          if (info/=0) return
          read(99, *) data1  !(1:N)
          read(99, *) data2  !(1:N)
      end if

      close(99)
      loadData2_CX = .true.
   end
 
!*****************************************************************
   logical function load2Data_CX(N, data1, data2, binMode, fileName)
      implicit none
      integer, intent(IN) :: N
      double complex, intent(OUT) :: data1(N), data2(N)
      logical, intent(IN) :: binMode
      character(LEN=*), intent(IN) ::  fileName

      integer :: i, info

      load2Data_CX = .false.
      if (binMode) then
          open(99,FILE=fileName, STATUS='OLD', FORM='UNFORMATTED',  &
               IOSTAT=info)
          if (info/=0) return
          do I = 1, N
             read(99) Data1(I), data2(I)
          end do
      else
          open(99, FILE=FILENAME,STATUS='OLD', IOSTAT=info)
          if (info/=0) return
          do I = 1, N
             read(99, *) data1(I), data2(I)
          end do
      end if
         
      close(99)
      load2Data_CX = .true.
    end 

!*************************************************************      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                     Print data                              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine printData_CX(N, DATA1)
      implicit none
      integer, intent(IN) :: N
      double complex,intent(IN) :: DATA1(N)

      print *
      print *, 'Print Complex Data. Data Size=', N

      print *, DATA1 !(1:N)
             
      end 
!**************************************************************
