!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!c          Diagonalize a Symmetric Matrix          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
! TVG: This function is a duplicate by name in mm.f90
! This appears to be redundant, so the name is changed to allow it to compile
! If enabling, remember to change the last lines to reflect the name
logical function DIAG_duplicate(JOBZ, N, H1, E0)
   implicit none
   character,intent(IN):: JOBZ
   integer, intent(IN) :: N
   double precision, intent(INOUT) :: H1(N,N)
   double precision, intent(OUT)   :: E0(N)
   
   integer, parameter :: dscale=3

   double precision   :: work(dscale*N) 
   integer :: lwork, info
   
   lwork = dscale*N

   call DSYEV(JOBZ,'U', N, H1, N, E0, WORK, LWORK, Info)

   diag_duplicate = .false.
   if (Info == 0)  diag_duplicate=.true.

   !ccccccccccccccccccccccccccccccccccccccccccccccccc
   !c  sometimes it is necessary to check the sign  c
   !c   of eigen-vector, for V and -V is the same   c
   !ccccccccccccccccccccccccccccccccccccccccccccccccc 

end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!c          Diagonalize a Symmetric Matrix          c
!cccccccccccccccccccccccccccccccccccccccccccccccccccc

logical function DIAGF(N, H1, E0, filename, srMode)
   implicit none
   integer, intent(IN) :: N
   double precision, intent(INOUT) :: H1(N,N)
   double precision, intent(OUT)   :: E0(N)
   character(len=*), intent(IN)    :: filename
   logical, intent(IN) :: srMode
   
   integer :: info
   logical :: diag

   if (srMode) then 
      diagF = diag('V', N, H1, E0)

      if (diagF) then
         open(99, FILE=filename, status='replace',  &
               form='unformatted', iostat=info)
         if (info == 0 ) then
             write(99) N, E0, H1
             close(99)
         else
             diagF = .FALSE.
         end if
      end if
   else
      open(99, FILE=filename, status='old',        &
               form='unformatted', iostat=info)
      if (info == 0 ) then
             write(99) info
             close(99)
             diagF = .TRUE.
         else
             diagF = .FALSE.
         end if
   end if
end 
