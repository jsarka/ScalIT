!
! Do VBR convergence testing 
!
program test_pso_lr2
   use presinc
   external fitVlr2
   
   print *, '***********************************'
   print *, '       1D DVR Calculation.         '
   print *, '***********************************'
  
   print *
   print *, ' Prepare data for further calculation' 
   print *, ' DVR preparation for BR'

   print *
   print *, ' Running ........'

   call run(fitVlr2)
   
   print *, '========= FINISHED =============='

end program



