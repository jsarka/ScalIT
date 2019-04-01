! have tested ADVRCx, ADVRCxF
program test_gausscx

   integer, parameter :: N=11
   double complex :: S1(N,N), X1(N,N)
   double precision :: X0(N)
   logical :: wrFlag=.true.
   integer :: opt
   character(len=*),parameter :: fname='./tgau.dat'

   double precision :: a, b, alpha

   a=1.0D0; b=2.0D0; alpha=70.0D0
   opt=3


   call GaussCX_XS(N,a,b,alpha,S1,X1)
   select case (opt)
   case (1)
      call ADVRCX(N,S1,X1,X0)

   case (2)
      print *, 'Save transformation matrices'
      wrFlag=.true.
      call ADVRCxF(N,S1,X1,X0,wrFlag,fname)
   
   case (3)
      wrFlag=.false.
      print *, 'Load transformation matrices'
      call ADVRCxF(N,S1,X1,X0,wrFlag,fname)
   end select

   print *, 'DVR points:'
   print *, X0

   print *, ' V^T*V-V*V^T:'
!   print *, MatMul(S1,X1)-MatMul(X1,S1)

   print *, ' V*VT:'
!   print *, MatMul(S1,X1)

end
