!
!
program test_jbech

   double precision :: lr(3),BR(3), cth(3)

   print *
   print *, ' Testing  program to exchange 2 particles in A3'
   print *, ' Input lr, Br, and cos(theta)'
   read(*,*) lr(1), BR(1), cth(1)

   lr(1)=ABS(lr(1))+1.0D-11; BR(1)=ABS(BR(1))+1.0D-11
   if (ABS(cth(1))>1.0D0) cth(1)=1.0D0

   call JBexch(BR(1),lr(1),cth(1), BR(2),lr(2),cth(2),BR(3),lr(3),cth(3));
   print *
   write(*,10) 'lr, A-BC:',lr

   call JBexch(BR(2),lr(2),cth(2), BR(3),lr(3),cth(3),BR(1),lr(1),cth(1));
   write(*,10) 'lr, B-CA:',lr

   call JBexch(BR(3),lr(3),cth(3), BR(1),lr(1),cth(1),BR(2),lr(2),cth(2));
   write(*,10) 'lr, C-AB:',lr

   call JBexch(BR(1),lr(1),cth(1), BR(2),lr(2),cth(2),BR(3),lr(3),cth(3));
   print *
   write(*,10) 'BR, A-BC:',Br

   call JBexch(BR(2),lr(2),cth(2), BR(3),lr(3),cth(3),BR(1),lr(1),cth(1));
   write(*,10) 'BR, B-CA:',Br

   call JBexch(BR(3),lr(3),cth(3), BR(1),lr(1),cth(1),BR(2),lr(2),cth(2));
   write(*,10) 'BR, C-AB:',BR

   call JBexch(BR(1),lr(1),cth(1), BR(2),lr(2),cth(2),BR(3),lr(3),cth(3));
   print *
   write(*,10)  'cos(Theta), A-BC:',cth

   call JBexch(BR(2),lr(2),cth(2), BR(3),lr(3),cth(3),BR(1),lr(1),cth(1));
   write(*,10)  'cos(Theta), B-CA:',cth

   call JBexch(BR(3),lr(3),cth(3), BR(1),lr(1),cth(1),BR(2),lr(2),cth(2));
   write(*,10)  'cos(Theta), C-AB:',cth
   print *, 'Finish the Testing!'
   print *

  10 format(A,  1x, 3(F15.9, 1X))
end
