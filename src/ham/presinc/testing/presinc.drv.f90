program test_presinc

   use presinc

   external :: fitV

   call run(fitV)

end

subroutine fitV(N, R, VR)
   integer, intent(IN) :: N
   double precision, intent(IN) :: R(N)
   double precision, intent(OUT):: VR(N)

   VR(1:N)=0.0D0

end 



