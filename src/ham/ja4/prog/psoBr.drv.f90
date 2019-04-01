!
! Do VBR convergence testing 
!
program test_pso_BR
   use presinc
   external fitVBR
   
   print *, '***********************************'
   print *, '       1D DVR Calculation.         '
   print *, '***********************************'
  
   print *
   print *, ' Prepare data for further calculation' 
   print *, ' DVR preparation for BR'

   print *
   print *, ' Running ........'

   call run(fitVBR)
   
   print *, '========= FINISHED =============='

end program



