!
! Do VBR convergence testing 
!
program test_pso_lr1
   use presinc
   external fitVlr1
   
   print *, '***********************************'
   print *, '       1D DVR Calculation.         '
   print *, '***********************************'
  
   print *
   print *, ' Prepare data for further calculation' 
   print *, ' DVR preparation for lr1'

   print *
   print *, ' Running ........'

   call run(fitVlr1)
   
   print *, '========= FINISHED =============='

end program



