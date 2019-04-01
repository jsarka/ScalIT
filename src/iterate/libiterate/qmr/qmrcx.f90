!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                            QMR Function                                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 integer function QMR_CX(nTYPE, QMRMAX, QMRTOL, N, B, HXCX, PXCX, X, RES)
      implicit none
      integer, intent(IN)           :: nTYPE, QMRMAX, N
      double precision, intent(IN)  :: QMRTOL
      double complex,intent(IN)     :: B(N)   
      external :: HXCX, PXCX
      double complex,intent(OUT)    :: X(N)  
      double precision,intent(OUT)  :: RES    

!cccccccccccccc Temporary vectors      
      double complex :: R(N), UX(N), UR(N) 
      double complex :: Q(N), S(N), V(N), Z(N)   

      double complex ::  DOT_CX   ! X^T * Y   

      double precision  ::  BNRM, ZNRM, VNRM, VNRMPO, ZVRAT
      double complex    :: QSDOT, VZDOT, BETA, ETA     
      double precision  :: GAMAMO, GAMA, THETAMO, THETA 

      integer  :: i, info                    
      double precision :: RES1, RES2, RES3, XNRM  
      double complex, allocatable :: HX0(:)

!ccccccccccccccccccccccccccccccccccc
      QMR_CX = -1
      if (ABS(nTYPE) == 1) then    ! conv. using |HX-B|/|HX| 
          allocate(HX0(N), STAT=info)
          if (info /= 0) return
      end if

!cccccccccccccccccccccccccccccccccccccccc
!c     Initialization: R = V = B        c
!cccccccccccccccccccccccccccccccccccccccc     

      ZVRAT = 1.0d0
      
      if (nType <= 0)  X(1:N) = 0.0D0
    
      R(1:N) = B(1:N);     V(1:N) = B(1:N)             
      BNRM   = DSQRT(abs(dot_product(B,B)))   
      VNRM   = BNRM

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    Z = P * B ~ H^-1 * B, the preconditioner         c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call PXCX(N, B, Z)     
      ZNRM   = DSQRT(abs(dot_product(Z, Z)))  
      GAMAMO = 1.0D0;       ETA    = -1.0D0

!ccccccccccccccccccccccccccccc      
!c     Main Loop             c
!ccccccccccccccccccccccccccccc

      QMR_CX = 0
      MAIN_LOOP: do i = 1, QMRMAX          

         if ( (VNRM == 0.0d0) .OR. (ZNRM == 0.0d0))  exit  

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c   BP addition to reduce number of vectors and matrix-vector products.  c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ZVRAT = ZVRAT * ZNRM/VNRM

!ccccccccccccccccccccccccccccccccccccccccccccc
!c     Normalize V and Z, compute VZDOT      c
!ccccccccccccccccccccccccccccccccccccccccccccc
         V = V / VNRM ;       Z = Z / ZNRM           
         VZDOT = DOT_CX(N, V, Z)  
         if (abs(VZDOT)==0.0D0)   exit 

!ccccccccccccccccccccccccccccccc
!c        Update Q             c
!ccccccccccccccccccccccccccccccc
         if ( i == 1) then
            Q(1:N) = Z(1:N)         
         else
            Q(1:N) = Z(1:N) - (VNRM*VZDOT)/QSDOT * Q(1:N)
         end if

!ccccccccccccccccccccccccccccccc         
!c    Update S = ZVRAT H Q     c
!ccccccccccccccccccccccccccccccc
         call HXCX(N, Q, S)
         S(1:N) = ZVRAT * S(1:N)    

!ccccccccccccccccccccccccccccccccccc
!c      Update QSDOT and BETA      c
!ccccccccccccccccccccccccccccccccccc
         QSDOT = DOT_CX(N, Q, S)  
         BETA  = QSDOT / VZDOT      

         if ((abs(QSDOT)==0.0D0) .OR. (abs(BETA)==0.0d0)) &   
            exit 

!ccccccccccccccccccccccccccccccc         
!c      Update V, VNRMPO       c
!ccccccccccccccccccccccccccccccc

         V(1:N) = S(1:N) - BETA * V(1:N)  
         VNRMPO = DSQRT(abs(dot_product(V,V)))   

!cccccccccccccccccccccccccccccccccccc
!c     Update Z = P^(-1) V/ZVRAT    c   
!cccccccccccccccccccccccccccccccccccc
         call PXCX(N, V, Z)              
         Z(1:N) = Z(1:N) / ZVRAT       

!ccccccccccccccccccccccccccccccccccccccccccc
!c    Update ZNRM, THETA, GAMA, ETA        c
!ccccccccccccccccccccccccccccccccccccccccccc
         
         ZNRM = DSQRT(abs(dot_product(Z, Z))) 
         THETA = VNRMPO/( GAMAMO * abs(BETA) )
         GAMA = 1.0d0/DSQRT(1.0d0+THETA*THETA)
         ETA = -ETA*VNRM*GAMA*GAMA/(BETA*GAMAMO*GAMAMO)

         if (GAMA == 0.0D0) exit    
!cccccccccccccccccccccccccccccccc         
!c      Update X and R          c
!cccccccccccccccccccccccccccccccc
         if (i == 1) then           
            UX(1:N) = ETA*ZVRAT * Q(1:N) 
            UR(1:N) = ETA * S(1:N)       
         else
            UX(1:N) = (ETA*ZVRAT)* Q(1:N) + ((THETAMO*GAMA)**2)*UX(1:N)
            UR(1:N) = ETA*S(1:N) + ((THETAMO*GAMA)**2)*UR(1:N)
         end if
                 
         X(1:N) = X(1:N) + UX(1:N)   
         R(1:N) = R(1:N) - UR(1:N)
   
!cccccccccccccccccccccccccccccccccccccccc 
!c     Relative residual processing     c
!cccccccccccccccccccccccccccccccccccccccc
         ! convergence criterion
         select case (ABS(nTYPE))     
         case (1)       ! |H*X-B|/|X| 
               call HXCX(N, X, HX0) 
               XNRM = DSQRT(abs(dot_product(HX0, HX0)))
               HX0(1:N) = HX0(1:N) - B(1:N)
               RES = DSQRT(abs(dot_product(HX0, HX0)))/XNRM

         case (2)       ! Full: |R|/|B|, |UX|/|X|, |UR|/|R|
               RES1 = DSQRT(abs(dot_product(R, R)))/BNRM 
               RES2 = DSQRT(abs(dot_product(UX,UX)/dot_product(X, X))) 
               RES3 = DSQRT(abs(dot_product(UR,UR)/dot_product(R, R))) 
               RES = max(RES1, RES2, RES3)

         case DEFAULT   ! Simple: |R|/|B|
               RES = DSQRT(abs(dot_product(R, R)))/BNRM  
                   ! RES = NORM(R,N)/BNRM, 
         end select

         if (i/200*200 ==i)    write(*,10) i, RES

         if (abs(RES) < QMRTOL) then
            QMR_CX = i;       exit
         else
            QMR_CX = -i  
         end if

!cccccccccccccccccccccccccccccccccccccccc
!c    Update variables on two levels    c
!cccccccccccccccccccccccccccccccccccccccc
         VNRM = VNRMPO;       GAMAMO = GAMA
         THETAMO = THETA

      end do  MAIN_LOOP

      if (allocated(HX0))    deallocate(HX0)

 10   format('      QMR Iter:',I8,3x,'RES:',2X,E20.15)     

 end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

