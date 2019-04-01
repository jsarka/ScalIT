!
! Do VBR convergence testing 
!
program test_pso_BR
   use presinc
   external fitVlr
   
   print *, '***********************************'
   print *, '       1D DVR Calculation.         '
   print *, '***********************************'
  
   print *
   print *, ' Prepare data for further calculation' 
   print *, ' DVR preparation for BR'

   print *
   print *, ' Running ........'

   call run(fitVlr)
   
   print *, '========= FINISHED =============='

end program



