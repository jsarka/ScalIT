!ccccccccccccccccccccccccccccccccccccccccc
!c   IO routines to calculate wave3d     c
!ccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   Subroutines to calculate wave-function        c
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calSaveWF()
      double precision :: ct1, ct2, ct3
      integer :: myChoice
   
      call CPU_Time(ct1)
      print *
      print *, ' Computing wave functions ......'   

      myChoice = kNum
      if ((kNum>=0).AND.(kNum<=mMax)) myChoice=0

      if (gType) then
          select case (myChoice)
          case (-2)
             call calWFP1()
          case (-1)
             call calWFS1()
          case (0)
             call calWFM1(kNum)
          case default
             call calWFS1()
          end select
      else
          select case (myChoice)
          case (-2)
             call calWFP2()
          case (-1)
             call calWFS2()
          case (0)
             call calWFM2(kNum)
          case default
             call calWFS2()
          end select
      end if

      call CPU_Time(ct2)
      print *, ' CPU Time to calculate wave function:', ct2-ct1

      print *
      print *, ' Save wave function in file:',fOut
      call saveWF()
      call CPU_Time(ct3)

      print *, ' CPU Time to save wave function:', ct3-ct2
   end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine saveWF()
   integer :: i, j, k

   if (sType) then
      open(99, file=fOut, form='unformatted', status='replace')
      write(99) NState, NSInd(1:NState),NR(1:3)
      write(99) myR1, myR2, myR3
      write(99) P0
      close(99)
   else
      open(99, file=fOut, form='formatted', status='replace')
      do i = 1, NR(1)
         do j = 1, NR(2)
            do k = 1, NR(3)
               write(99,10)  myR1(i),myR2(j),myR3(k),p0(1:NState,i,j,k), p0(1:NState,i,j,k)**2
            end do
         end do
      end do
      close(99)
   end if

  10 format(3(F10.6),1x,20(E15.9,2x))

end subroutine 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!cccccccccccccccccccccccccccccccccccccccc
!c          Read input parameters       c
!cccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine readW3StdIO()
    
   call myreadW3(5)

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine readW3File(filename)
   character(len=*), intent(IN) :: filename

   open(99, File=Filename, status='old')

   call myreadW3(99)    

   close(99)

end subroutine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
